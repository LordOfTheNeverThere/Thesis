from __future__ import annotations #  postpone evaluation of annotations
import pandas as pd
from itertools import combinations
import numpy as np
import math
from sklearn.metrics import auc
from scipy.stats import pearsonr
import pickle
import gzip
from resources import *




class MatrixData:
    def __init__(self, filepath: str = None, data: pd.DataFrame = None, **readerKwargs):
        self.data = data
        assert filepath or (
            data is not None), 'There should be either a filepath or data'

        if filepath:
            self.data: pd.DataFrame = pd.read_csv(filepath, **readerKwargs)
            
        elif data is not None:
            self.data: pd.DataFrame = data.copy()
        
    def __str__(self) -> str:
        return str(self.data)

    def write(self, filepath:str):

        with gzip.open(filepath, 'wb') as f:
            pickle.dump(self,f)
        f.close()

    def query(self, query: str, inplace : bool = False) -> pd.DataFrame|None:


        if not inplace:
            return self.data.query(query).copy()
        else:
            self.data = self.data.query(query).copy()

    def compare(self, other: MatrixData, querySelf: str, queryOther: str, key: str = 'PPI') -> pd.DataFrame:
        """Query two PairwiseCorrMatrices with independent queries and get the merged dataframe as result

        Args:
            other (PairwiseCorrMatrix): another PairwiseCorrMatrix object
            querySelf (str): Query to apply on the self object
            queryOther (str): Query to apply on the other object
            key (str, optional): In what column name should the merge be happening on, acting as consesual index. Defaults to 'PPI'.

        Returns:
            pd.DataFrame: _description_
        """

        left: pd.DataFrame = self.query(querySelf).copy()
        right: pd.DataFrame = other.query(queryOther).copy()

        return left.merge(right, on=key)
class ppiDataset(MatrixData):

    def __init__(self, filepath:str = None, data: pd.DataFrame = None, proteinLabels: list = [], **readerKwargs):


        super().__init__(filepath, data, **readerKwargs)
        self.proteinLabels = proteinLabels
        self.ppis = set()


    def getPPIs(self, dataset:str) -> set:
        """Get the curated and observed ppis of a certains external PPI dataset

        Args:
            dataset (str): name of the dataset. It can either be 'corum', 'biogrid' or 'string'

        Returns:
            set: Set of tuples with proteinA and proteinB as putative PPI pairs
        """

        data = self.data.copy()
        ppiSet = None
        allowedDatasets = ['corum', 'biogrid', 'string']

        assert dataset in allowedDatasets, f"dataset not supported use one of the following 'corum', 'biogrid', 'string', got: {dataset}"

        if dataset == allowedDatasets[0]:

            def combinationsOfProteins(complx):

                if ';' in list(complx['subunits(Gene name)']):
                    complx['proteinTuple'] = list(combinations(
                        complx['subunits(Gene name)'].split(';'), 2))

                    return complx

            data = data.apply(
                lambda complx: combinationsOfProteins(complx), axis=1)
            ppiList = list(data.dropna()['proteinTuple'])
            ppiSet = {item for sublist in ppiList for item in sublist}

        elif (dataset == allowedDatasets[1]):

            # Filter Biogrid for certain parameters

            # Only allow ppis documented by physical interaction      
            # Filter out Homedymers which are not object of our study
            data = data.query("`Experimental System Type` == 'physical' and `Official Symbol Interactor A` != `Official Symbol Interactor B`").copy()

            data['proteinTuple'] = list(zip(data['Official Symbol Interactor A'], data['Official Symbol Interactor B']))
            ppiSet = set(data['proteinTuple'])

            

        elif (dataset == allowedDatasets[2]):
            data['proteinTuple'] = list(zip(data['proteinA'], data['proteinB']))
            ppiSet = set(data['proteinTuple'])

        
        self.ppis = ppiSet
        return ppiSet


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
    
    def tlsResidues(self) -> ResiduesMatrix:

        proteomics = self.data.copy()
        tlsResList = []
        

        for columnX in proteomics:
            for columnY in proteomics:

                if (columnX == columnY): #We dont want to calculate the residues of homo Pairwise 'Correlations'
                    continue

                index= columnX + ';' + columnY
                X = proteomics.loc[:,columnX].dropna(axis=0) #Get X and Y data
                Y = proteomics.loc[:,columnY].dropna(axis=0)
                samplesInCommon = X.index.intersection(Y.index)

                if len(samplesInCommon) < 5: # We have no important information from protein protein interactions with less than 5 coocorence
                    continue

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
                residues = pd.DataFrame(residues,columns=[index])
                tlsResList.append(residues)


        tlsResData = pd.concat(tlsResList, join='outer', sort=False)
        print(tlsResData)

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


class PairwiseCorrMatrix(MatrixData):

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, ** readerKwargs):

        super().__init__(filepath, data, **readerKwargs)

        self.corrCumSum = None
        self.indexes = None
        self.auc = None
        self.label = None

    def __str__(self):

        print = super().__str__()
        print = print + '\n' +str(self.auc) + '\n' + str(self.label)

        return print

    def addGroundTruth(self, ppis: set, externalDatasetName: str):
        """Append the binary values of a putative PPI, from an external dataset (e.g Corum), to our pairwise correlation Dataframe

        Args:
            ppis (ppiDataset): ppiDataset of the external ppi dataset used
            data (pd.DataFrame): Pairwise correlation dataframe
            externalDatasetName (str): Name to give to the binary column holding the truth value of an PPI is seen in that external Dataset
            filepath (str): The name to give to the final file, with the added binary column

        Returns:
            _type_: Data with added column
        """
        data = self.data.copy()
        def addExternalTrueY(model):


            found = 0
            [proteinA, proteinB] = model.name.split(';')
            # In my implementation the ppis have (A,B) but not (B,A), they are combinations
            ppiAB: tuple = (proteinA, proteinB)
            ppiBA: tuple = (proteinB, proteinA)

            if ppiAB in ppis or ppiBA in ppis:
                found = 1

            model[externalDatasetName] = found

            return model[externalDatasetName]

        data[externalDatasetName] = data.apply(
            axis=1, func=lambda model: addExternalTrueY(model))


        self.data = data

        return data
    
    def aucCalculator(self, yColumnName:str, label:str):
        """Adds the value of AUC of the Recall curve using a specified external PPI dataset with yColumnName

        Args:
            yColumnName (str): Name of the df column where there is the truth value of the existence or not of the PPI in the reported external PPI dataset
            label (str): Text which will show up as label next to the value of the AUC, e.g 'Baseline Auc == 0.9' 
        """
        pairwiseCorr = self.data 

        self.corrCumSum = np.cumsum(
            pairwiseCorr[yColumnName]) / np.sum(pairwiseCorr[yColumnName])
        
        self.indexes = np.array(pairwiseCorr.reset_index().index) / \
            pairwiseCorr.shape[0]
        self.auc = auc(self.indexes, self.corrCumSum)

        if not label:
            self.label = f"(AUC {self.auc:.2f})"
        
        self.label = label + f" (AUC {self.auc:.2f})"

    
class DrugResponseMatrix(MatrixData):
    """Class interface and methods for the drug response data"""

    def __init__(self, filepath: str=None, data: pd.DataFrame=None, **readerKwargs):
        super().__init__(filepath, data, **readerKwargs)


    def binrise(self, deathThresh: int = 3, inplace: bool = True):
        """Creates a Binary representation of the Drug response matrix, where 0 means no efficacy and 1 efficacy. This is done bye calculating a threshold, 50% of the natural log of the [max screen] 


        Args:
            deathThresh (int, optional): Number required of 'responsive' cell lines for a certain drug to be considered in the final set. Defaults to 3.
        """
        data = self.data.copy()
        maxScreenConc = pd.read_csv(PATH + '/externalDatasets/drugMaxScreenConcentration.csv', index_col='Unnamed: 0')
        maxScreenConc.index.name = 'drug'
        
        data = data.merge(maxScreenConc, on='drug') # We only work with drugs for which we have a max screen concentration and IC50
        # With 50% of the natural log of the max Screen Concentration we can create an empirically proven threshold, 
        # where generelly IC50 values greater than that mean that The drug has no efficacy, and below drug has efficancy, 
        # or in other words the cell line doesn't or does 'respond to a drug'
        data['efficacyThreshold']= data['MAX_CONC'].apply(lambda row: math.log(row * 0.5)) 
        data.drop(columns=['MAX_CONC'], inplace=True)


        # Create mask that see if all columns are less than the values in the threshold col, and convert them to ints so the bool's become either 0 or 1
        data:pd.DataFrame = data.apply(lambda row: row < row['efficacyThreshold'], axis=1).astype(int)

        relevantDrugs = (data.sum(axis=1) >= deathThresh) # Condition that only accounts for drugs that kill at leats 3 cell lines

        if inplace:
            self.data = data.loc[relevantDrugs].drop(columns=['efficacyThreshold']) # apply condition (We lose 130 drugs, 697 -> 567) 26/4/23
        else:
            return data.loc[relevantDrugs].drop(columns=['efficacyThreshold'])
        
class ResiduesMatrix(MatrixData):

    def __init__(self, filepath: str=None, data: pd.DataFrame=None, **readerKwargs):
        super().__init__(filepath, data, **readerKwargs)

    



        




        

