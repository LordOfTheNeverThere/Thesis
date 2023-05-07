import numpy as np
import pandas as pd
from scipy.special import stdtr
import matplotlib.pyplot as plt
import time as t
from resources import *

def getGLSCorr(proteinData: ProteinsMatrix, pValues: bool = True, listCovMatrix:list[pd.DataFrame] = None, coefColumnName :str = 'glsCoefficient') -> PairwiseCorrMatrix:

    proteinData = proteinData.data.copy()

    if listCovMatrix is not None:
        # The matrix used for the covariance is not the one used as an X in the linear regression, 
        # we used this to allow the use of the genomic matrix as the matrix where the cov of each sample would be calculated, 
        # since the tecidual bias of proteomic expression would be more present in the genome, 
        # since this layer is upstream and more connected to the information the cell has from its diferentiation, 
        # so it becomes easier to distinguish samples by tissue, by using its genomic data. 
        # Therefore calculating the cov of genomics.csv could give a more correct value on the true covariation of two samples 
        # and correct the possible correlation of PPIs that are simply due to the diferent baseline expression of proteins belonging in the ppi

        listCovMatrix.append(proteinData)
        samplesInCommon = pd.concat(listCovMatrix, axis=1, join='inner').index #Conactenate all matrices and select their indices which would be the common ones
        listCovMatrix.pop(-1) # Remove recently added proteinData to the list

        proteinData = proteinData.loc[samplesInCommon, :] 
        proteinDataMean = proteinData.fillna(proteinData.mean())
        proteinDataMean.dropna(axis=1, inplace=True) #We delete all columns with nan values because after all this preprocessing they must be completetly empty columns
        covMatrix = np.zeros((len(samplesInCommon), len(samplesInCommon)))

        for dataForCov in listCovMatrix:
            
            dataForCov = dataForCov.loc[samplesInCommon, :]
            dataForCov.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #How we are handling missing data, there should be at leats 20% of missingness for a collumn to be dropable
            dataForCov = dataForCov.fillna(dataForCov.mean())
            print(dataForCov)
            # calculate covariance matrix in order to see the covariace between samples, and notice tecidual patterns codified in the samples
            dataForCov = np.cov(dataForCov)
            print(dataForCov)
            covMatrix = covMatrix + dataForCov
            print(covMatrix)



    else: # The matrix used for the covariance is the same as that used for X in linear regression, we are whitening while taking into account the covariation of our data in X
        
        proteinData.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column
        proteinDataMean = proteinData.fillna(proteinData.mean())
        dataForCov = proteinDataMean
        # calculate covariance matrix in order to see the covariace between samples, and notice tecidual patterns codified in the samples
        covMatrix = np.cov(dataForCov)

    proteinNames = proteinDataMean.columns.str.split(' ').str.get(0).to_numpy()
    proteinNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]


    #invert it the covariance matrix
    covMatrix = np.linalg.inv(covMatrix)

    # Decompose it with Cholesky, returning the lower triangular matrix of the positive definite matrix covMatrix, because cov(x1,x2) == cov(x2,x1)
    cholsigmainvMean = np.linalg.cholesky(covMatrix)


    # Whittening transformation, we codify our data into a space where each the variance of each covariate is the same and equal to one, 
    # so we are kind like normalising it, in fact that's exactly what we are doing ~ N(0,I) As they call it warping...
    warpedProteinsMean = proteinDataMean.T.values @ cholsigmainvMean

    # The Mean of the variances of a sample across all gene symbols is used as a beta0 for linear regression
    warpedIntereceptMean = cholsigmainvMean.T.sum(axis=0)


    def linear_regression(warped_screens, warped_intercept):
        GLS_coef = np.empty((len(warped_screens), len(warped_screens)))
        GLS_se = np.empty((len(warped_screens), len(warped_screens)))
        ys = warped_screens.T
        for proteinIndex in range(len(warped_screens)):
            
            X = np.stack((warped_intercept, warped_screens[proteinIndex]), axis=1)
            if np.any(np.isnan(X)) :
                print(proteinIndex)
                print('\n')
                print(X)
            coef, residues = np.linalg.lstsq(X, ys, rcond=None)[:2]
            df = warped_screens.shape[1] - 2
            GLS_coef[proteinIndex] = coef[1]
            GLS_se[proteinIndex] = \
                np.sqrt(np.linalg.pinv(X.T @ X)[1, 1] * residues / df)
        return GLS_coef, GLS_se

    GLS_coef, GLS_se = linear_regression(warpedProteinsMean, warpedIntereceptMean)

    df = warpedProteinsMean.shape[1] - 2

    #   Construct new PairwiseGLSCoefs Matrix


    glsCoefs = GLS_coef[np.triu_indices(GLS_coef.shape[0], k=1)]
    if pValues: #We might not want to add pValues so that we use less memory
        GLS_p = 2 * stdtr(df, -np.abs(GLS_coef / GLS_se))
        np.fill_diagonal(GLS_p, 1)
        glsPValues = GLS_p[np.triu_indices(GLS_p.shape[0], k=1)]
        pairwiseCorrData = pd.DataFrame(
            {coefColumnName: glsCoefs, 'p-value': glsPValues}, index=proteinNames)
    else:
        pairwiseCorrData = pd.DataFrame(
            {coefColumnName: glsCoefs}, index=proteinNames)
        
    pairwiseCorrData.sort_values(by='glsCoefficient', ascending=False, inplace=True) #Sorting the pairwise corr by the higest beta coeficient as proxy to PPI
    
    pairwiseCorrData = PairwiseCorrMatrix(None, pairwiseCorrData)
    pairwiseCorrData.data.index.name = 'PPI'

    return pairwiseCorrData


if __name__ == '__main__':
    gls: PairwiseCorrMatrix = utils.read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    glsRNA: PairwiseCorrMatrix = utils.read(PATH + '/datasetsTese/glsPairwiseCorrTranscriptCov.pickle.gz')
    glsRNAProt: PairwiseCorrMatrix = utils.read(PATH + '/datasetsTese/glsPairwiseCorrTranscriptProteinCov.pickle.gz')
    corum: ppiDataset = utils.read(PATH + '/externalDatasets/corum.pickle.gz')

    print(gls)
    print(glsRNA)
    print(glsRNAProt)

    gls.addGroundTruth(corum.ppis, 'corum')
    glsRNA.addGroundTruth(corum.ppis, 'corum')
    glsRNAProt.addGroundTruth(corum.ppis, 'corum')
    gls.aucCalculator('corum', 'gls Model')
    glsRNA.aucCalculator('corum', 'gls covRNA')
    glsRNAProt.aucCalculator('corum', 'gls covRNAProtein')

    gls.write(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    glsRNA.write(PATH + '/datasetsTese/glsPairwiseCorrTranscriptCov.pickle.gz')
    glsRNAProt.write(PATH + '/datasetsTese/glsPairwiseCorrTranscriptProteinCov.pickle.gz')

    print(gls)
    print(glsRNA)
    print(glsRNAProt)


    
    

