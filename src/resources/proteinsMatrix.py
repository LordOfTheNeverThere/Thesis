import numpy as np
import pandas as pd
from matrixData import MatrixData
from scipy.special import stdtr
from scipy.stats import pearsonr



class ProteinsMatrix(MatrixData):

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, **readerKwargs):

        super().__init__(filepath, data, **readerKwargs)

    def pearsonCorrelations(self, columnName: str, counting: bool = True, pValue: bool = True) -> PairwiseCorrMatrix:
        """Calculate the pearson correlations and corresponding p-value, displaying them in a pairwise manner, returning an instance of the PairwiseCorrMatrix class

        Args:
            fileName (str): Name of the file where the datastructure will be stored
            columnName (str): Name given to the df column with the correlation metric
            counting (bool, optional): Should we count the number of times the two proteins appear in the same sample?. Defaults to True.
            pValue (bool, optional): Should we add the p-value, level of statistical significane of our correlation to the final data structure. Defaults to True.

        Returns:
            PairwiseCorrMatrix: final data structure with all the information regarding pairwise correlations
        """

        data = self.data.copy()
        # Get list with the names of every PPI
        proteinNames = data.columns.str.split(' ').str.get(0).to_numpy()
        ppiNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]
        # Correlation Matrix
        pearsonCorrMatrix = data.corr(method='pearson', numeric_only=True)

        pairwiseCorrData = pearsonCorrMatrix.to_numpy()[np.triu_indices(pearsonCorrMatrix.shape[0], k=1)]

        pairwiseCorrData = pd.DataFrame({columnName: pairwiseCorrData}, index=ppiNames)
        pairwiseCorrData.index.names=['PPI']

        if counting:

            # Co-occorance Matrix
            # Get 1 and 0 depending on if there is a value or not in a specific spot
            coOccuranceMatrix = (data/data).fillna(0).astype(int)
            # Simple linear algebra to get the co-occurance values
            coOccuranceMatrix = coOccuranceMatrix.T.dot(coOccuranceMatrix)
            coOccuranceMatrix = coOccuranceMatrix.to_numpy()[np.triu_indices(coOccuranceMatrix.shape[0], k=1)]
            coOccuranceData = pd.DataFrame({'counts': coOccuranceMatrix}, index=ppiNames)
            coOccuranceData.index.names = ['PPI']

            pairwiseCorrData = pairwiseCorrData.merge(
                coOccuranceData, on='PPI', how='left')
            
        if pValue:

            def pearsonPValues(data:pd.DataFrame = None)-> pd.DataFrame|None:

                pValuesMatrix = data.corr(method=lambda x, y: pearsonr(x, y)[1])
                pairwisePValues =pValuesMatrix.to_numpy()[np.triu_indices(pValuesMatrix.shape[0], k=1)]
                pairwisePValues = pd.DataFrame({'pValue': pairwisePValues}, index=ppiNames)
                pairwisePValues.index.names = ['PPI']

                return pairwisePValues['pValue']
            
            pairwiseCorrData['pValue'] = pearsonPValues(pearsonCorrMatrix)
            

        pairwiseCorrData.sort_values(
            by=columnName, ascending=False, inplace=True)


        return PairwiseCorrMatrix(None,pairwiseCorrData.dropna()) #There will be NAN correlations between proteins which do not appear simultaneously in at least two cell lines
    
    def tlsResidues(self, pairwiseCorr:PairwiseCorrMatrix) -> ResiduesMatrix:

        pairwiseCorr = pairwiseCorr.data
        proteomics = self.data.copy()
        tlsResList = []
        ppis =  pairwiseCorr.loc[pairwiseCorr['corum'] == 1].loc[pairwiseCorr['glsCoefficient'] > 0.8].index 
        #Get the ppis that are most likely true ppis, so that we can analyse what samples do not correspond to the correlation, 
        # are farthest from the linear regession line, hence, have greatest TLS and so are samples of interest where the PPI likely. 
        # So it would be interesting to se afterwards if that sample has a responsiveness to a drug all the other samples do not meaning 
        # we are in a presence of a PPI that might be correlated to a feature, a certain drug responsiveness
        

        for ppi in ppis:
            proteinA = ppi.split(';')[0]
            proteinB = ppi.split(';')[1]

            X = proteomics.loc[:,proteinA].dropna(axis=0) #Get X and Y data
            Y = proteomics.loc[:,proteinB].dropna(axis=0)
            samplesInCommon = X.index.intersection(Y.index)

            if len(samplesInCommon) < 5: # We have no important information from protein protein interactions with less than 5 coocorence
                continue
                print(f"{ppi} was not used for calculating the tls Redidues Matrix because it did not have at least 5 samples")

            X=X.loc[samplesInCommon] #Locate samples that both have some value (not nan)
            Y=Y.loc[samplesInCommon]

            meanX = X.mean()
            meanY = Y.mean()

            meanErrorX = X - meanX #Start Calculating quantities to minimise for the tlsCoef Calculation
            meanErrorY = Y - meanY

            meanSqErrorX = meanErrorX ** 2
            meanSqErrorY = meanErrorY ** 2

            u = meanSqErrorX.sum()
            v = meanSqErrorY.sum()
            r = (meanErrorX * meanErrorY).sum()
            w = v - u

            tlsCoef = (w + (w**2 + r**2)**0.5) #Calculating tls Coefficient
            tlsCoef = tlsCoef/r

            intercept = meanY - (tlsCoef * meanX) #Intercept of linear fit
            predY = intercept + (tlsCoef * X)
            residues = Y - predY # TLS Residues in absolute val
            residues = pd.DataFrame(residues,columns=[ppi])
            tlsResList.append(residues)


        tlsResData = pd.concat(tlsResList, join='outer', sort=False)
        return ResiduesMatrix(None,tlsResData)

    def getGLSCorr(self, pValues: bool = True, listCovMatrix:list[pd.DataFrame] = None, coefColumnName :str = 'glsCoefficient') -> PairwiseCorrMatrix:
        """Get the GLS coeficents between each Protein X and Y, where X != Y, these will measure the correlation between each protein. 
        But this method differs from the pearsonCorrelations since it has underneath a GLM where the covariance matrix can be any specified.
        This covariance matrix will transform both X and y of the proteinData.data. By default this covariance matrix is calculated with proteinData.data
        Where we get the covariance between samples, as a similarity measure between samples. This tweak is speacilly important if our residuals
        correlate with X, meaning we are not in the normal OLS case.

        Args:
            pValues (bool, optional): Add the pValues of each gls Coefficient to the output data. Defaults to True.
            listCovMatrix (list[pd.DataFrame], optional):List of matrices to use to calculate covariance, if only one is required insert [matrix]. Defaults to None.
            coefColumnName (str, optional): Name to appear on the Dataframe' Column of glsCoefs. Defaults to 'glsCoefficient'.

        Returns:
            PairwiseCorrMatrix: Data structure with all above information
        """
        proteinData = self.data.copy()

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

