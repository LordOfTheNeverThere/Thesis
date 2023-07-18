from __future__ import annotations #  postpone evaluation of annotations
import pandas as pd
from itertools import combinations, product, permutations, repeat, chain
import numpy as np
import math
from sklearn.metrics import auc
import pickle
import gzip
from scipy.special import stdtr
from scipy.stats import pearsonr
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.diagnostic import het_breuschpagan, het_white
from statsmodels.stats.multitest import multipletests
from statsmodels.multivariate.pca import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import normalize, QuantileTransformer, StandardScaler
from scipy.stats import chi2, shapiro
import seaborn as sns
import matplotlib.pyplot as plt
import time as t
from typing import Iterable,Any
import multiprocessing as mp
from resources import *


def ols(Y,X):
    """Give the OLS regression results for a given set of features and a response variable, like R Summary

    Args:
        Y (_type_): Data with dependent variables
        X (_type_): Data with independent variables

    Returns:
        OLS: Regressor object, used for predictions for example
        RegressionResultsWrapper:  Results of the regression, like betas and all types of test statistics
    """
    X = sm.add_constant(X)
    regressor = sm.OLS(Y, X)
    results = regressor.fit()
    print(results.summary())

    return regressor, results

def covMatrixAnalysis(data:pd.DataFrame)-> tuple[float, float]:
    """Gives a summary of the covariance matrix in terms of the % of  and the mean of the non diagonal elements,
    useful to see if the observations are independent or not, if the mean of the non diagonal elements is close to zero and 
    that of the diagonal is close to 1 then the observations are independent

    Args:
        data (pd.DataFrame): Dataframe from which the covariance matrix will be calculated

    Returns:
        tuple[float, float]: (Mean of the diagonal, mean of the non diagonal elements of the covariance matrix)
    """
    cov = np.cov(data)
    upperTriangular = np.triu(cov, k=1)
    upperTriangular = upperTriangular.flatten()
    upperTriangular = [round(x, 8)==0 for x in upperTriangular]
    nonDiagonalMetric:float = ((sum(upperTriangular))/len(upperTriangular)) * 100
    #I only did for the upper triangular because the covariance matrix is symmetric
    diagonal = np.diag(cov)
    diagonalMetric = [round(x, 4) == 1 for x in diagonal]

    diagonalMetric = ((sum(diagonalMetric))/len(diagonalMetric)) * 100  

    return diagonalMetric, nonDiagonalMetric

def calcR(XY) -> tuple[float, float]:
    X,Y = XY
    corr, pValue = pearsonr(X, Y)
    return corr, pValue


def calcMahalanobis(y:pd.DataFrame, data: pd.DataFrame, cov:pd.DataFrame=None):

    y_mu = (y - np.mean(data, axis=0)).T # In Liner Algebra The covariates are usually row vectors in dataframes these are usually column vectors
    if not cov:
        cov = np.cov(data.values.T)
    inv_covmat = np.linalg.inv(cov)
    left = np.dot(y_mu.T, inv_covmat)
    mahal = np.dot(left, y_mu).diagonal()
    pValue = 1 - chi2.cdf(mahal, 3)

    return np.sqrt(mahal), pValue

class MatrixData:
    def __init__(self, filepath: str = None, data: pd.DataFrame = None, **readerKwargs):
        self.data = data
        self.filepath = filepath

        if filepath is not None:
            self.data: pd.DataFrame = pd.read_csv(filepath, **readerKwargs)
            
        elif data is not None:
            self.data: pd.DataFrame = data.copy()
        
    def __str__(self) -> str:
        return str(self.data)

    def write(self, filepath:str = None):

        if filepath is not None:
            self.filepath = filepath

        if self.filepath is not None:
            filepath = self.filepath

        assert filepath is not None, 'No filepath provided'

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
    
    def getUniqueSetValues(self, feature: str):
        """Returns a set of unique values from a feature of a dataframe

        Args:
            feature (str): The column name to extract the unique set of values

        Returns:
            set: The uniqye set of values in a column of a Dataframe
            dict: Dictionary with the keys are the unique values in a column and the values(of the dict) as the number of
            occurances in of each value(of the feature)
        """

    
        data = self.data.copy()
        setOfValues = data[feature].unique()
        setOfValues = set(setOfValues)
        # occurancesDict = data.groupby(feature).count().to_dict()[feature]

        return setOfValues, None

class ppiDataset(MatrixData):

    def __init__(self, name:str, filepath:str = None, data: pd.DataFrame = None, **readerKwargs):


        super().__init__(filepath, data, **readerKwargs)
        self.name = name
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
                    complx['proteinTuple'] = list(permutations(
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

            ppis = list(zip(data['Official Symbol Interactor A'], data['Official Symbol Interactor B'])) + list(zip(data['Official Symbol Interactor B'], data['Official Symbol Interactor A']))
            ppiSet = set(ppis)

            

        elif (dataset == allowedDatasets[2]):
            ppis = list(zip(data['proteinA'], data['proteinB'])) + list(zip(data['proteinB'], data['proteinA']))
            ppiSet = set(ppis)

        
        self.ppis = ppiSet
        return ppiSet


class ProteinsMatrix(MatrixData):

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, **readerKwargs):

        super().__init__(filepath, data, **readerKwargs)

    def __str__(self) -> str:
        string =  super().__str__()

        try:

            for summary in self.normSummary:    
                string +=  "\n" + f"For a Global p-value of {summary[0]:.2f} and a threshold of {summary[1]} samples, {summary[2]:.2f}% of the proteins are not normally distributed"
            
            for summary in self.homoskeSummary:    
                string +=  "\n" + f"For a Global p-value of {summary[0]:.2f} and a threshold of {summary[1]} samples, {summary[2]:.2f}% in {summary[3]} PPIs of after linear regression don't follow the homoskedasticity assumption"

            return string
        
        except:
            return string

    def pearsonCorrelations(self, columnName: str, proteomicsType:str,thresholdInteraction:int = 5) -> PairwiseCorrMatrix:
        """Calculate the pearson correlations and corresponding p-value, displaying them in a pairwise manner, returning an instance of the PairwiseCorrMatrix class

        Args:
            columnName (str): Name given to the df column with the correlation metric
            proteomicsType (str): The type of proteomics data to use to create the pearson correlation matrix
            thresholdInteraction (int, optional): The minimum number of coincident samples for it to be considered a putative PPI

        Returns:
            PairwiseCorrMatrix: final data structure with all the information regarding pairwise correlations
        """
        data = self.data.copy()
        permutationsColumns = list(combinations(data.columns, 2))
        pairwiseCorr:dict = {'PPI':[],'proteinA': [],'proteinB': [], columnName: [], 'p-value':[], 'counts':[]}
        
        for (proteinAName, proteinBName) in permutationsColumns:

            proteinA = data[proteinAName].dropna(axis=0)
            proteinB = data[proteinBName].dropna(axis=0)
            samples  = proteinA.index.intersection(proteinB.index)
            proteinA = proteinA.loc[samples]
            proteinB = proteinB.loc[samples]

            count = len(proteinA)
            if count  < thresholdInteraction :
                continue
            
            
            (corr, pValue) = pearsonr(proteinA, proteinB)
            pairwiseCorr['PPI'].append(proteinAName + ';' + proteinBName)
            pairwiseCorr['proteinA'].append(proteinAName)
            pairwiseCorr['proteinB'].append(proteinBName)
            pairwiseCorr['pearsonR'].append(corr)
            pairwiseCorr['p-value'].append(pValue)
            pairwiseCorr['counts'].append(count)

        index = pairwiseCorr.pop('PPI')
        return PairwiseCorrMatrix(proteomicsType, None, pd.DataFrame(pairwiseCorr, index=index), ['pearsonR', 'p-value'], [False, True])
      
    def calculateResidues(self, ppis: Iterable[set[str]]) -> ResiduesMatrix:
    

        proteomics = self.data.copy()
        tlsResList = []
        correlationsTLSMahal = []

        # are farthest from the linear regession line, hence, have greatest TLS and so are samples of interest where the PPI likely is having some biomolecular role. 
        # So it would be interesting to see afterwards if that sample has a responsiveness to a drug all the other samples do not meaning 
        # we are in a presence of a PPI that might be correlated to a feature, a certain drug responsiveness

        for index, ppi in enumerate(ppis):
 
            proteinA = ppi.split(';')[0]
            proteinB = ppi.split(';')[1]

            X = proteomics.loc[:,proteinA].dropna(axis=0) #Get X and Y data
            Y = proteomics.loc[:,proteinB].dropna(axis=0)
            samplesInCommon = X.index.intersection(Y.index)

            if len(samplesInCommon) < 5: # We have no important information from protein protein interactions with less than 5 coocorence
                
                print(f"{ppi} was not used for calculating the tls Redidues Matrix because it did not have at least 5 samples")
                continue

            X=X.loc[samplesInCommon] #Locate samples that both have some value (not nan)
            Y=Y.loc[samplesInCommon]
            # TLS Residues
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
            residues = abs(Y - predY) # TLS Residues in absolute val

            # Malahanobis Distance
            proteinExpression = pd.concat([X,Y], axis=1)
            mahalDist, mahalPValues = calcMahalanobis(proteinExpression, proteinExpression, None)

            dfData = {(ppi, 'TLS'): residues,
                      (ppi,'malahanobis'): mahalDist,
                      (ppi,'mahalPValue'): mahalPValues,
                     }
            correlationsTLSMahal.append(pearsonr(residues, mahalDist)[0])
            residues = pd.DataFrame(dfData)
            if index == 0:
                tlsResData = residues
            else:
                tlsResList.append(residues)

        print('Statistical Description of the Pearson Correlation between TLS and Mahalanobis distance \n' + str(pd.DataFrame(correlationsTLSMahal).describe())) 
        tlsResData = tlsResData.join(tlsResList, how='outer')
        return ResiduesMatrix(None,tlsResData)
    

    @classmethod
    def whitening(cls, proteinData:ProteinsMatrix, covMatrix:ProteinsMatrix, saveIndexes:bool= False) -> tuple[Any, Any]:
        """Whitten the proteinData, so that each covariate has the same variance, which is equal to one

        Args:
            covMatrix (_type_):The Protein Data to calculate the covariance Matrix in order to withen the data

        Returns:
            tuple[Any, Any]: the warped X to be used in the linear regression and the intercept, which is the mean of the variances of a sample across all gene symbols
        """        
        #invert it the covariance matrix
        proteinData = proteinData.data
        covMatrix = covMatrix.data
        samplesCommon = proteinData.index.intersection(covMatrix.index)
        if len(samplesCommon) < len(proteinData.index):
            print(f"We have lost {len(proteinData.index) - len(samplesCommon)} samples ")

        proteinData = proteinData.loc[samplesCommon]
        covMatrix = covMatrix.loc[samplesCommon]
        
        covMatrix = np.cov(covMatrix)
        covMatrix = np.linalg.inv(covMatrix)

        # Decompose it with Cholesky, returning the lower triangular matrix of the positive definite matrix covMatrix, because cov(x1,x2) == cov(x2,x1)
        cholsigmainvMean = np.linalg.cholesky(covMatrix)


        # Whittening transformation, we codify our data into a space where each the variance of each covariate is the same and equal to one, 
        # so we are kind like normalising it, in fact that's exactly what we are doing ~ N(0,I) As they call it warping...
        warpedProteinsMean = proteinData.T.values @ cholsigmainvMean

        # The intercept is the sum of the choleski decomposition matrix, when the sum equals to one that sample is independent from all the others
        warpedIntereceptMean = cholsigmainvMean.T.sum(axis=0)

        if saveIndexes:
            warpedProteinsMean = pd.DataFrame(warpedProteinsMean.T, columns=proteinData.columns, index=proteinData.index)



        return warpedProteinsMean, warpedIntereceptMean

    def getGLSCorr(self, proteomicsType:str, pValues: bool = True, listCovMatrix:list[pd.DataFrame] = None, coefColumnName :str = 'glsCoefficient') -> PairwiseCorrMatrix:
        """Get the GLS coeficents between each Protein X and Y, where X != Y, these will measure the correlation between each protein. 
        But this method differs from the pearsonCorrelations since it has underneath a GLM where the covariance matrix can be any specified.
        This covariance matrix will transform both X and y of the proteinData.data. By default this covariance matrix is calculated with proteinData.data
        Where we get the covariance between samples, as a similarity measure between samples. This tweak is speacilly important if our residuals
        correlate with X, meaning we are not in the normal OLS case.

        Args:
            proteomicsType (str): The type of proteomics data to used to calculate the GLS Coefficients
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
                # calculate covariance matrix in order to see the covariace between samples, and notice tecidual patterns codified in the samples
                dataForCov = np.cov(dataForCov)
                covMatrix = covMatrix + dataForCov




        else: # The matrix used for the covariance is the same as that used for X in linear regression, we are whitening while taking into account the covariation of our data in X
            
            proteinData.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column
            proteinDataMean = proteinData.fillna(proteinData.mean())
            dataForCov = proteinDataMean
            # calculate covariance matrix in order to see the covariace between samples, and notice tecidual patterns codified in the samples
            covMatrix = np.cov(dataForCov)


        proteinNames = proteinDataMean.columns.str.split(' ').str.get(0).to_numpy()
        proteinNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]
        proteinsA, proteinsB = [proteinPair.split(";")[0] for proteinPair in proteinNames], [proteinPair.split(";")[1] for proteinPair in proteinNames]

        warpedProteinsMean, warpedIntereceptMean = ProteinsMatrix.whitening(self, self)


        def linear_regression(warped_screens, warped_intercept):
            GLS_coef = np.empty((len(warped_screens), len(warped_screens)))
            GLS_se = np.empty((len(warped_screens), len(warped_screens)))
            ys = warped_screens.T
            for proteinIndex in range(len(warped_screens)):
                
                X = np.stack((warped_intercept, warped_screens[proteinIndex]), axis=1)
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
                {'proteinA':proteinsA, 'proteinB': proteinsB, coefColumnName: glsCoefs, 'p-value': glsPValues}, index=proteinNames)
        else:
            pairwiseCorrData = pd.DataFrame(
                {'proteinA':proteinsA, 'proteinB': proteinsB, coefColumnName: glsCoefs}, index=proteinNames)
        
        pairwiseCorrData = PairwiseCorrMatrix(proteomicsType,None, pairwiseCorrData, [coefColumnName, 'p-value'], [False, True])
        pairwiseCorrData.data.index.name = 'PPI'

        return pairwiseCorrData
    
    def getGLSR(self, proteomicsType:str, coefColumnName :str = "GLM's R") -> PairwiseCorrMatrix:

        proteinData = self.data.copy()
        # Get protein names in pair and individually
        proteinNames = proteinData.columns.str.split(' ').str.get(0).to_numpy()
        proteinNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]
        proteinsA, proteinsB = [proteinPair.split(";")[0] for proteinPair in proteinNames], [proteinPair.split(";")[1] for proteinPair in proteinNames]
        #Transform data into Cholesky Space with Whitening
        withenedProteins, withenedInterecept = ProteinsMatrix.whitening(self, self, saveIndexes=True)
        #Get X and Y for Pearson correlation, now with the withened data
        X = [withenedProteins[proteinA].to_numpy() for proteinA in proteinsA]
        Y = [withenedProteins[proteinB].to_numpy() for proteinB in proteinsB]
        

        # Calculate Pearson correlation in parallel
        results = {
            'PPI': proteinNames,
            'proteinA': proteinsA, 
            'proteinB':proteinsB,
            coefColumnName: [], 
            'p-value': []}


        with mp.Pool(CPUS) as process:
            pararelResults = process.map(calcR, list(zip(X, Y)))

        for corr, pValue in pararelResults:

            results[coefColumnName].append(corr)
            results['p-value'].append(pValue)

        # Extract the results
        index = results.pop('PPI')
        pairwiseCorrData = pd.DataFrame(results, index=index)
        pairwiseCorr = PairwiseCorrMatrix(proteomicsType, None, pairwiseCorrData, [coefColumnName, 'p-value'], [False, True])

        return pairwiseCorr


        

    def shapiroWilksTest(self, thresh: int = 5, globalPVal:float = 0.01) -> None:

        """Performs the shappiro wilks test for each protein present and stores it in self.normTest and in self.normSummary

        Args:
            thresh (int): Minimum number of samples to perform the test
            globalPVal (float): Global p-value threshold to reject the null hypothesis

        Returns:
            Nonetype: None
        """

        data = self.data.copy()
        
        shapiroResults = {}


        for protein in data:

            proteinData = data[protein].dropna()    

            if len(proteinData) >= thresh:
                stat, pVal = shapiro(proteinData)
                
            
            else:
                stat, pVal = np.nan, np.nan

            shapiroResults[protein] = {'stat': stat, 'p-value': pVal}
        
        shapiroResults = pd.DataFrame(shapiroResults).T
        self.normTest = shapiroResults
        shapiroResults = shapiroResults.dropna(axis=0)
        pValues = shapiroResults['p-value']
        rejected, correctedPVals, _, _ = multipletests(pValues, alpha=globalPVal, method='fdr_bh')
        numNonNormal = np.sum(rejected)
        ratioNonNormal = (numNonNormal/shapiroResults.shape[0]) * 100 # The smaller the pValue the more likely it is that the data is not normally distributed, we thence reject the null hypothesis that the data is normally distributed

        try: # If the atribute doesn't exist we create it
            self.normSummary.add((globalPVal,  thresh, ratioNonNormal))
        except:
            self.normSummary = set()
            self.normSummary.add((globalPVal,  thresh, ratioNonNormal))
        
        print(self.normSummary)

    def whiteTest(self,  thresh: int = 5, ppis:set = None , globalPVal:float = 0.05) -> None:
        """Executes  the White test (where H_0 is the homoskedasticity of the residuals, thence if the residuals are invariant to the change of x, 
        meaning that for y~x, x explains most of the variability, leaving no confounding factor out of the regression equation), for a global pValue, and with PPIs with at leats thresh samples

        Args:
            thresh (int, optional): The minimum number of samples required to calculate the test. Defaults to 5.
            globalPVal (float, optional): p-value subject to the bonferroni correction. Defaults to 0.01.

        Returns:
            Nonetype: None
        """

        data = self.data.copy()
        if ppis is None:
            ppis: permutations[tuple[str, str]] = permutations(data.columns, 2)

        whiteResults = {}

        for x,y in ppis:


            #Getting and Processing Data
            xData = data[x].dropna()
            yData = data[y].dropna()
            samplesCommon = xData.index.intersection(yData.index)

            if len(samplesCommon) >= thresh:
                
                xData = xData.loc[samplesCommon]
                yData = yData.loc[samplesCommon]
                
                #Fitting the model (Classical Linear regression)
                regressor = LinearRegression()
                regressor.fit(xData.values.reshape(-1,1), yData.values.reshape(-1,1)) 
                yPred = regressor.predict(xData.values.reshape(-1,1))
                #Calculating the residuals
                residuals = yData.values.reshape(-1,1) - yPred
                #Add intercept to the data
                xData = pd.DataFrame({'x':xData, 'intercept':np.ones(len(xData))})

                #Calculating the White test
                stat, pValue,_,_ = het_breuschpagan(residuals, xData)

                whiteResults[(x,y)] = {'stat': stat, 'p-value': pValue}
        
        whiteResults = pd.DataFrame(whiteResults).T.reset_index(names=['proteinA', 'proteinB']) # This allows for a compatible merging with PaiwiseCorrMatrix objects
        self.homoskeTest = whiteResults
        numPPIs = whiteResults.shape[0]
        whiteResults = whiteResults.dropna(axis=0)
        pValues = whiteResults['p-value']
        rejected, _, _, _ = multipletests(pValues, alpha=globalPVal, method='fdr_bh')
        numHeteroske = np.sum(rejected)
    
        ratioHeteroske = (numHeteroske/numPPIs) * 100 # The smaller the pValue the more likely it is that the residuals are heteroskedastic, we thence reject the null hypothesis that the residuals are invariant when regressed with x, homoskedastic
        try: # If the atribute doesn't exist we create it
            self.homoskeSummary.add((globalPVal,  thresh, ratioHeteroske))
        except:
            self.homoskeSummary = set()
            self.homoskeSummary.add((globalPVal,  thresh, ratioHeteroske, numPPIs))

        print(self.homoskeSummary)


    def independenceSamples(self) -> None:

        data = self.data.copy()

        diagonalMean, nonDiagonalMean = covMatrixAnalysis(data)
        whitenedData, _ = ProteinsMatrix.whitening(self, self, True)
        whiteDiagonalMean, whiteNonDiagonalMean = covMatrixAnalysis(whitenedData)

        self.samplesIndep = f"Percentage of 1's in Diagonal: {diagonalMean}\nPercentage of 0's in non diagonal: {nonDiagonalMean}\nPercentage of 1's in Diagonal, after whitening: {whiteDiagonalMean}\nPercentage of 0's in non diagonal, after whitening: {whiteNonDiagonalMean}"
        print(self.samplesIndep)

    def plotPxPyDrug(self, drug:str, ppi:str, drugResponse: DrugResponseMatrix, filepath:str, **annotationArgs):


        drugResponse = drugResponse.binrise(inplace=False).T # The drug response matrix is binarised
        samplesCommon = self.data.index.intersection(drugResponse.index) # We only want to plot the samples that are in both matrices
        assert len(samplesCommon) > 0, 'There are no samples in common between the protein data and the drug response data'
        drugResponse = drugResponse.loc[samplesCommon, drug]
        proteinData = self.data.loc[samplesCommon, :]

        if len(ppi.split('-')) > 0:
            ppi = ppi.split('-')[0]

        pxName = ppi.split(';')[0]
        px = proteinData[pxName]
        pyName = ppi.split(';')[1]
        py = proteinData[pyName]

        data = pd.DataFrame({'px': px, 'py': py, 'drugResponse': drugResponse})
        colors = {0: 'green', 1: 'red'}
        #change dot size onm the scatter plot
        plt.figure(figsize=(15, 15))	
        plt.scatter(data['px'], data['py'], c=data['drugResponse'].map(colors), label=[colors.values(), colors.keys()])
        plt.title('Protein expression \n with Drug Response, >50% [drugMax]')
        plt.xlabel(str(pxName))
        plt.ylabel(str(pyName))
        if annotationArgs is not None:
            plt.annotate(**annotationArgs)
        legend = [plt.Line2D([0], [0], marker='o', color=colors[key], label='Drug Response = ' + str(key)) for key in colors]
        plt.legend(handles = legend, fontsize=20, framealpha=0.2)
        plt.savefig(filepath)

        plt.close()

    def plotPxPyDrugContinous(self, drug:str, ppi:str, drugResponse: DrugResponseMatrix, filepath:str, **annotationArgs):
        """Scatter Plot with the protein expression of two proteins and the drug response of a drug, in a continous manner, unlike plotPxPyDrug.
        Additionally, the plot can be annotated with the arguments passed to the function. And the Drug Response will be represented with a colorbar and
        and an overlaying Kernel Distribution Estimation.

        Args:
            drug (str): name of the drug
            ppi (str): Protein-Protein Interaction (Px;Py)
            drugResponse (DrugResponseMatrix): Drug Response Matrix
            filepath (str): filepath to save the plot
        """
        drugResponse = drugResponse.data.T # The drug response matrix is binarised
        samplesCommon = self.data.index.intersection(drugResponse.index) # We only want to plot the samples that are in both matrices
        assert len(samplesCommon) > 0, 'There are no samples in common between the protein data and the drug response data'
        drugResponse = drugResponse.loc[samplesCommon, drug]
        proteinData = self.data.loc[samplesCommon, :]

        if len(ppi.split('-')) > 0:
            ppi = ppi.split('-')[0]

        pxName = ppi.split(';')[0]
        pyName = ppi.split(';')[1]
        plottingData = proteinData[[pxName, pyName]]
        plottingData = plottingData.join(drugResponse, how = 'inner')
        #standardize the data 
        plottingData = pd.DataFrame(StandardScaler().fit_transform(plottingData) ,columns=plottingData.columns, index=plottingData.index)


        plt.figure(figsize=(10, 10))
        scatter = sns.scatterplot(data=plottingData, x=pxName, y=pyName, hue=drug, palette="viridis", alpha=1, edgecolor='none', s=15)
        # Add Colour Map
        sm = plt.cm.ScalarMappable(cmap="viridis")
        sm.set_array([])
        scatter.get_legend().remove()
        scatter.figure.colorbar(sm, label='Drug Response')

        plt.title('Protein expression \n with Drug Response')
        plt.xlabel(str(pxName))
        plt.ylabel(str(pyName))

        if annotationArgs is not None:
            plt.annotate(**annotationArgs)

        plt.savefig(filepath)
        plt.close()



    
    def plotPxPySample(self, ppis: list[str], filepath: str, sampleProp: str) -> None:
        """Plots the px against py proteins for each ppi in ppis, and colors the samples according to the sampleProp
        """

        data = self.data.copy()
        samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
        samplesCommon = data.index.intersection(samplesheet.index)
        if len(samplesCommon) < data.shape[0]:
            print('We have lost samples for proteomics since the samplesheet did not completely overlap with our samples, check for an updated samplesheet')
        data = data.loc[samplesCommon, :]
        samplesheet = samplesheet.loc[samplesCommon, :].fillna('NaN')
        uniquePropValues = list(MatrixData(None, samplesheet).getUniqueSetValues(sampleProp)[0])
        # if len(uniquePropValues) > 400:
        #     print(
        #         f'There are more than 400 unique values for {sampleProp}, which is not very informative, please choose a different sampleProp')
        #     return None
        colors = {value: color for value, color in zip(uniquePropValues, sns.color_palette('dark', len(uniquePropValues)))}
        numPPIs = len(ppis)

        # Setting up plot axes

        fig, ax = plt.subplots(numPPIs, 1, figsize=(15, 15*numPPIs))

        for index, ppi in enumerate(ppis):

            pxpyNames = ppi.split(';')
            pxName, pyName = pxpyNames[0], pxpyNames[1]
            px, py = data[pxName], data[pyName]
            samplesCommon = px.index.intersection(py.index)
            # Get only relevant samples and data
            px, py = px.loc[samplesCommon], py.loc[samplesCommon]
            
            try:
                samplePropData =  samplesheet.loc[samplesCommon, sampleProp]
            except:
                print("No Samples in Commmon between pxpy and feature")
                continue
            

            plottingData = pd.DataFrame(
                {'px': px, 'py': py, sampleProp: samplePropData})

            ax[index].scatter(plottingData['px'], plottingData['py'], c=plottingData[sampleProp].map(colors), label=[colors.values(), colors.keys()])
            ax[index].set_title('test')
            ax[index].set_xlabel(pxName)
            ax[index].set_ylabel(pyName)
            legend = [plt.Line2D([0], [0], marker='.', color=colors[key], label=f'{sampleProp} = ' + str(key)) for key in colors]
            ax[index].legend(handles=legend, fontsize=8, framealpha=0.2)

        # Save
        fig.savefig(filepath)
        plt.close()

    
    def PCA(self, filepath:str='', numPC:int = 10, gls:bool = False, factorsName:str='', **pcaKwargs) -> PCA:
        """Generate the PCA of the protein data

        Args:
            filepath (str, optional): File name and directory where 
            the scree and PC cumulative variance plot will be stored. Defaults to ''.
            numPC (int, optional): number of principal components used. Defaults to 10.
            gls (bool, optional): If the PCA class will use the gls parameter 
            from the documentation it states that "Flag indicating to implement a two-step GLS 
            estimator where in the first step principal components are used to estimate residuals, 
            and then the inverse residual variance is used as a set of weights to estimate the final 
            principal components. Setting gls to True requires ncomp to be less then the min of the 
            number of rows or columns". Defaults to False.

        Returns:
            PCA: The returned object will be the PCA object with all the relevant atributes
        """        
        # Fit PCA
        pca = PCA(self.data, ncomp=numPC, gls=gls, missing='drop-row', **pcaKwargs)

        # Extract explained variance of each PC
        explainedVar = pca.eigenvals / np.sum(pca.eigenvals)

        # Calculate cumulative explained variance
        cumulativeVar = np.cumsum(explainedVar)
        
        # Plot Scree plot along with cumulative explained variance
        scree = sns.barplot(x=np.arange(1, numPC + 1), y=explainedVar, color='blue') # Scree plot, individual explained variance of PCA
        sns.lineplot(x=np.arange(0, numPC), y=cumulativeVar, color='red', ax=scree) # Cumulative explained variance of PCA
        # Change labels of axes
        scree.set(xlabel='Principal Component', ylabel='Explained Variance')
        #show figure
        plt.show()
        # Save figure
        if filepath != '':
            plt.savefig(filepath, bbox_inches='tight')
        if factorsName != '':
            #Change the name of the columns of the factors and scores given by the PCA Object
            newColumns = [f'{factorsName}{i}' for i in range(1, numPC + 1)]
            pca.factors.columns = newColumns
            pca.scores.columns = newColumns

        return pca
        





class PairwiseCorrMatrix(MatrixData):

    def __init__(self, proteomicsType:str, filepath: str = None, data: pd.DataFrame = None, proxies:list[str] = ['p-value', 'coef'], ascendings:list[bool] = [True, False],** readerKwargs):
        """Initialize the PairwiseCorrMatrix object, where we calculate the correlation between all the proteins in the data, and their p-value
        Additionally we can use one of these values to calculate the auc of that proxy metric towards recalling PPIs in a specific yColumn Dataset, the ground truth

        Args:
            filepath (str, optional): filepath where Datafrane is stored to instatiate the object. Defaults to None.
            data (pd.DataFrame, optional):Dataframe to instatiate Object. Defaults to None.
            aucs (_type_, optional): dict with name of columns to which we calculate the auc and value. Defaults to None.
        """
        super().__init__(filepath, data, **readerKwargs)
        self.yColumn = []
        self.proxies= proxies
        self.ascendings = ascendings
        self.labels = {proxy:{} for proxy in self.proxies}
        self.corrCumSums = {proxy:{} for proxy in self.proxies}
        self.indexes = {proxy:{} for proxy in self.proxies}
        self.aucs = {proxy:{} for proxy in self.proxies} 
        self.proteomicsType = proteomicsType




    def __str__(self):

        print = super().__str__()
        for proxy, proxyDict in self.aucs.items():

            print = print + '\n' + str(proxy) + '\n'

            for yColumn, aucVal in proxyDict.items():

                print = print + str(self.labels[proxy][yColumn])
        

        return print
    

    
    def addGroundTruth(self, ppis: set, externalDatasetName: str):
        """Append the binary values of a putative PPI, from an external dataset (e.g Corum), to our pairwise correlation Dataframe

        Args:
            ppis (ppiDataset): ppiDataset of the external ppi dataset used
            data (pd.DataFrame): Pairwise correlation dataframe
            externalDatasetName (str): Name to give to the binary column holding the truth value of an PPI is seen in that external Dataset


        Returns:
            _type_: Data with added column
        """
        data = self.data.copy()

        data[externalDatasetName] = [int((pA,pB) in ppis) for pA,pB in zip(data['proteinA'], data['proteinB'])]


        self.data = data
        self.yColumn.append(externalDatasetName)

        return data

    @classmethod
    def addGroundTruths(cls, insts:Iterable[PairwiseCorrMatrix]):
        from .utils import read

        #get all ppis
        corum:ppiDataset = read(PATH + '/external/ppiDataset/corum.pickle.gz')
        biogrid:ppiDataset = read(PATH + '/external/ppiDataset/biogrid.pickle.gz')
        stringLow:ppiDataset = read(PATH + '/external/ppiDataset/string150.pickle.gz')
        stringMedium:ppiDataset = read(PATH + '/external/ppiDataset/string400.pickle.gz')
        stringHigh:ppiDataset = read(PATH + '/external/ppiDataset/string700.pickle.gz')
        stringHighest:ppiDataset = read(PATH + '/external/ppiDataset/string900.pickle.gz')

        ppis = [corum, biogrid, stringLow, stringMedium, stringHigh, stringHighest]

        for inst in insts:

            

            print(inst.proteomicsType)
            if 'proteinA' not in set(inst.data.columns) or 'proteinB' not in set(inst.data.columns):
                inst.data['proteinA'] = [ppi.split(';')[0] for ppi in inst.data.index]
                inst.data['proteinB'] = [ppi.split(';')[1] for ppi in inst.data.index]
                

            for ppi in ppis:
                #if the ppi set was already added to the dataframe, skip it
                if ppi.name in inst.yColumn:
                    continue
                inst.addGroundTruth(ppi.ppis, ppi.name)
            
            #Save what was done!
            inst.write()

    
    def aucCalculator(self, yColumnName:str, proteomicsType:str, proxyColumn:str, ascending:bool):
        """Adds the value of AUC of the Recall curve using a specified external PPI dataset with yColumnName

        Args:
            yColumnName (str): Name of the df column where there is the truth value of the existence or not of the PPI in the reported external PPI dataset
            label (str): Text which will show up as label next to the value of the AUC, e.g 'Baseline Auc == 0.9' 
            proxyColumn (str): Name of the column of the statistical meausure to quantify the probability of PPI existence
            ascending(bool): should the proxyColumn be ordered ascending or not for best AUC calculation
        """
        pairwiseCorr = self.data 

        
        pairwiseCorr.sort_values(by=proxyColumn, ascending=ascending, inplace=True) # We sort rows by the smallest to greatest pValues
        self.corrCumSums[proxyColumn][yColumnName] = np.cumsum(pairwiseCorr[yColumnName]) / np.sum(pairwiseCorr[yColumnName])
        self.indexes[proxyColumn][yColumnName] = np.array(pairwiseCorr.reset_index().index) / pairwiseCorr.shape[0]
        self.aucs[proxyColumn][yColumnName] = auc(self.indexes[proxyColumn][yColumnName], self.corrCumSums[proxyColumn][yColumnName]) # update aucs dict to have a new auc for a specific proxy column

        # if not label: #if the user did not insert any label default it
        #     self.labels[proxyColumn] = f"(AUC {proxyColumn} {self.aucs[proxyColumn]:.2f})"
        
        self.labels[proxyColumn][yColumnName] =  f" ({proteomicsType} proteomics using {proxyColumn} ⇒ AUC:{self.aucs[proxyColumn][yColumnName]:.2f} recalling {yColumnName})"

    def aucsCalculator(self, proteomicsType:str, proxyColumnList:list[str], ascendingList:list[bool], filepath:str = None):
        """Get the auc's of the recall curve for each proxy column (coef, p-value, etc) and for each external PPI dataset (corum, biogrid, etc)

        Args:
            proteomicsType (str): Name of the PairwiseCorrMatrix object
            proxyColumnList (list[str]): List of proxies for PPI existence
            ascendingList (list[bool]): List of booleans to indicate if the proxyColumn should be ordered ascending or not for best AUC calculation
            filepath (str, optional): Filepath where to save the pairwiseCorrMatrix Object. Defaults to None.
        """        

        for aucIndex in range(len(proxyColumnList)):
            for yColumn in self.yColumn:
                self.aucCalculator(yColumn, proteomicsType, proxyColumnList[aucIndex], ascendingList[aucIndex])
    
        if filepath is not None:
            self.write(filepath)
        print(self)

    @classmethod    
    def heatmap(cls, insts:list[PairwiseCorrMatrix], columns: list[str], bins:int, proteomics:ProteinsMatrix, filepath:str, title:str):
        """Creates a heatmap of the pairwise correlation between two columns of two PairwiseCorrMatrix objects, in order to see the order of missing values in a specifc range of both Pairwse Proxy correlation range

        Args:
            insts (list[PairwiseCorrMatrix]): Two PairwiseCorrMatrix objects
            columns (list[str]): Names of the proxy columns for each PairwiseCorrMatrix object
            bins (int):     Number of bins to use for the heatmap
            proteomics (ProteinsMatrix): Proteomics object to get the msv values
            filepath (str): Path to save the heatmap
            title (str): Title of the heatmap

        Returns:
            _type_: Dataframe with the values of the heatmap, and the number of PPIs in each bin
        """                
        dfs = [instance.data[column].copy() for instance, column in zip(insts, columns)]
        df = pd.concat(dfs, join='inner', axis=1)

        # Bin the two series that make up the dataframe with equal bins and return the intervals of each bin used
        df['bin0'] = pd.qcut(df[columns[0]], bins, precision=2)
        df['bin1'] = pd.qcut(df[columns[1]], bins, precision=2)
        intervals0 = sorted(df['bin0'].unique())
        intervals1 = sorted(df['bin1'].unique())

        heatmapData = pd.DataFrame()
        heatmapNumPPIs = pd.DataFrame()

        for interval0, interval1 in list(product(intervals0, intervals1)):


            colData = df.loc[df['bin0'] == interval0]
            rowData = df.loc[df['bin1'] == interval1]
            ppisCommon = set(colData.index.intersection(rowData.index)) # What are the ppis in common by the two queries
            mvs = 0 #missing values counter
 
            for ppi in ppisCommon: # count missing values
                
                proteinA = ppi.split(';')[0]
                proteinB = ppi.split(';')[1]
                mv =  proteomics.data[[proteinA, proteinB]].isna().sum().sum()  
                mvs =  mvs + mv
                numPPIs = len(ppisCommon)

                if numPPIs == 0:
                    mvsPerPPI = 0
                else:
                    mvsPerPPI = mvs / numPPIs #Standardise Mv in a query by the total number of ppis belonging to that query

            heatmapData.loc[str(interval0),str(interval1)] = mvsPerPPI
            heatmapNumPPIs.loc[str(interval0),str(interval1)] = numPPIs

        plt.figure(figsize=(8,8))
        sns.heatmap(heatmapData, annot=True, cmap='YlOrRd', fmt=".1f")
        plt.xlabel('Pearson R')
        plt.xticks(rotation=0)
        plt.ylabel('$β_{GLS}$')
        plt.title(title)

        plt.savefig(filepath)

        plt.close()

        #Number of PPIS per range heatmap
        plt.figure(figsize=(8,8))
        sns.heatmap(heatmapNumPPIs, annot=True, cmap='YlOrRd', fmt=".0f")
        plt.xlabel('Pearson R')
        plt.xticks(rotation=0)
        plt.ylabel('$β_{GLS}$')
        plt.title(title + '\n' +' Number of PPIs')

        plt.savefig(filepath.split('.')[0] + '#PPIs' + filepath.split('.')[1])
        
        return heatmapData, heatmapNumPPIs            

    @classmethod
    def getAucs(cls,instances:Iterable[PairwiseCorrMatrix]):
        """Calculates the Aucs for a list of PairwiseCorrMatrix objects

        Args:
            instances (Iterable[PairwiseCorrMatrix]): List of PairwiseCorrMatrix objects
        """ 
        for instance in instances:
            instance.aucsCalculator(instance.proteomicsType, instance.proxies, instance.ascendings, instance.filepath)


    
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
        maxScreenConc = pd.read_csv(PATH + '/external/drugMaxScreenConcentration.csv', index_col='Unnamed: 0')
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

    def getLinearModel(self, drugResponse: DrugResponseMatrix, samplesheet:pd.DataFrame, residualsType:str='TLS')->ResidualsLinearModel:

        X = self.data.copy()
        Y: pd.DataFrame = drugResponse.data.copy().T # Samples should be rows and not columns
        Y = Y.fillna(Y.mean())

        
        confoundingFactors = samplesheet[['tissue']].dropna(axis=0, how='any')
        confoundingFactors['lung'] = (confoundingFactors['tissue'] == 'Lung').astype(int)
        confoundingFactors['hymCellLine'] = (confoundingFactors['tissue'] == 'Haematopoietic and Lymphoid').astype(int)
        confoundingFactors = confoundingFactors.drop(columns=['tissue'])

        regressor = ResidualsLinearModel(Y, X, confoundingFactors, residualsType=residualsType)
        regressor.fit_matrix()

        return regressor



class GeneralLinearModel(MatrixData):
    """Authored by professor Emanuel Gonçalves, this call allows fo the computation of the linear regression coefs for a set of features. 
    Do check if the response variable is Normally distributed, since this is nor a GLM.
    There is no logging and covariate normalisation should be done prior to creating the object.
    """
    
    

    def __init__(
        self,
        Y,
        X,
        M,
        M2=None,
        fit_intercept=True,
        copy_X=True,
        n_jobs=4,
        verbose=1,
        addNoise:bool = False,
        filepath: str=None, 
        data: pd.DataFrame=None, 
        **readerKwargs
    ):
        super().__init__(filepath, data, **readerKwargs)  
        self.samples = set.intersection(
            set(Y.index),
            set(X.index),
            set(M.index),
            set(Y.index) if M2 is None else set(M2.index),
        )
        self.samples = list(self.samples) # make sure it's a list, because indexing by sets is deprecated
        self.X = X.loc[self.samples]
        self.X = self.X.loc[:, self.X.count() > (M.shape[1] + (1 if M2 is None else 2))]
        self.X = pd.DataFrame(StandardScaler().fit_transform(self.X), index=self.X.index, columns=self.X.columns)
        self.X_ma = np.ma.masked_invalid(self.X.values)
        

        self.Y = Y.loc[self.samples]
        self.Y = self.Y.loc[:, self.Y.std() > 0]


        self.M = M.loc[self.samples]

        self.M2 = M2.loc[self.samples, self.X.columns] if M2 is not None else M2

        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.n_jobs = n_jobs

        self.verbose = verbose
        self.addNoise = addNoise
        # self.log = logging.getLogger("Crispy")

    def model_regressor(self):
        regressor = LinearRegression(
            fit_intercept=self.fit_intercept,
            copy_X=self.copy_X,
            n_jobs=self.n_jobs,
        )
        return regressor

    @staticmethod
    def loglike(y_true, y_pred):
        nobs = len(y_true)
        nobs2 = nobs / 2.0

        ssr = np.power(y_true - y_pred, 2).sum()

        llf = -nobs2 * np.log(2 * np.pi) - nobs2 * np.log(ssr / nobs) - nobs2

        return llf

    @staticmethod
    def multipletests_per(
        associations, method="fdr_bh", field="pval", fdr_field="fdr", index_cols=None
    ):
        index_cols = ["y_id"] if index_cols is None else index_cols

        d_unique = {tuple(i) for i in associations[index_cols].values}

        df = associations.set_index(index_cols)

        df = pd.concat(
            [
                df.loc[i]
                .assign(fdr=multipletests(df.loc[i, field], method=method)[1])
                .rename(columns={"fdr": fdr_field})
                for i in d_unique
            ]
        ).reset_index()

        return df

    def fit_matrix(self):
        lms = []

        for x_idx, x_var in enumerate(self.X):
            # if self.verbose > 0:
            #     self.log.info(f"LM={x_var} ({x_idx})")


            # Mask NaNs
            x_ma = np.ma.mask_rowcols(self.X_ma[:, [x_idx]], axis=0)

            # Build matrices
            x = self.X.iloc[~x_ma.mask.any(axis=1), [x_idx]]
            y = self.Y.iloc[~x_ma.mask.any(axis=1), :]


            # Covariate matrix (remove invariable features and add noise)
            m = self.M.iloc[~x_ma.mask.any(axis=1), :]
            if self.M2 is not None:
                m2 = self.M2.iloc[~x_ma.mask.any(axis=1), [x_idx]]
                m = pd.concat([m2, m], axis=1)
            # m = m.loc[:, m.std(numeric_only = True) > 0]
            if self.addNoise:
                m += np.random.normal(0, 1e-6, m.shape)


            # Fit covariate model
            lm_small = self.model_regressor().fit(m, y)
            lm_small_ll = self.loglike(y, lm_small.predict(m))

            # Fit full model: covariates + feature
            lm_full_x = np.concatenate([m, x], axis=1)
            lm_full = self.model_regressor().fit(lm_full_x, y)
            betasFeature =  lm_full.coef_[:,-1]
            meanBetasCovariates = lm_full.coef_[:, 0:-1]

            lm_full_ll = self.loglike(y, lm_full.predict(lm_full_x))

            # Log-ratio test
            lr = 2 * (lm_full_ll - lm_small_ll)
            lr_pval = chi2(1).sf(lr)
            rSquared = lm_full.score(lm_full_x, y)
            #if Rsquared is significant print the association that is significant
            if rSquared > 0.65:
                print(f"Association: {x_var} and {y.columns} is significant with R^2 = {rSquared}")

            # Assemble + append results
            res = pd.DataFrame(
                dict(
                    y_id=y.columns,
                    x_id=x_var,
                    n=y.attrs["nan_mask"].loc[y.columns, x.index].sum(1) if "nan_mask" in y.attrs else len(x),
                    beta=betasFeature,
                    lr=lr.values,
                    covs=m.shape[1],
                    pval=lr_pval,
                    fdr=multipletests(lr_pval, method="fdr_bh")[1],
                    rSquared = rSquared
                )
            )
            
            smallerBetas = pd.DataFrame(meanBetasCovariates)
            smallerBetas.columns = [str(col) + 'Beta' for col in m.columns]
            res = pd.concat([res, smallerBetas], axis=1)
            lms.append(res)

        lms = pd.concat(lms, ignore_index=True).sort_values("pval")
        self.data = lms #Regression Results, of each X towwards every Y

        return lms

    @staticmethod
    def lm_residuals(y, x, fit_intercept=True, add_intercept=False):
        # Prepare input matrices
        ys = y.dropna()

        xs = x.loc[ys.index].dropna()
        xs = xs.loc[:, xs.std() > 0]

        ys = ys.loc[xs.index]

        if ys.shape[0] <= xs.shape[1]:
            return None

        # Linear regression models
        lm = LinearRegression(fit_intercept=fit_intercept).fit(xs, ys)

        # Calculate residuals
        residuals = ys - lm.predict(xs) - lm.intercept_

        # Add intercept
        if add_intercept:
            residuals += lm.intercept_

        return residuals
    
class TLSRegression():
    """Implements Total Least Squares Regression"""
    def __init__(self, Y:pd.DataFrame, X:pd.DataFrame, copy_X = True, fitIntercept=False, standardise:bool = True):

        self.Y = Y
        self.X = X
        self.samples = set.intersection(
            set(Y.index),
            set(X.index),
        )
        self.samples = list(self.samples) # make sure it's a list, because indexing by sets is deprecated
        self.X = X.loc[self.samples]
        if standardise:
            self.X = pd.DataFrame(StandardScaler().fit_transform(self.X), columns=self.X.columns, index=self.X.index) # Standardize X
        self.X_ma = np.ma.masked_invalid(self.X.values)

        self.Y = Y.loc[self.samples]

        self.fitIntercept = fitIntercept

        
    
    @staticmethod
    def tlsRegression(X:pd.Series|pd.DataFrame, Y:pd.Series|pd.DataFrame, fitIntercept=True):
        """Calculates the TLS regression of X on Y.

        Args:
            X (_type_): Covariates, excluding intercept
            Y (_type_): Response variable

        Returns:
           residuals, betas, predY, predX _type_: 
           (The regression's Residuals calculated with Forbenious norm of the errors of both X and Y, regression coefficients (includes beta0), predicted Y values, predicted X Values)
        """

        features = ['intercept'] + list(X.columns)
        samples = list(X.index)
        assert samples == list(Y.index), "X and Y must have the same samples"

        X = X.to_numpy() # Convert to numpy array
        Y = Y.to_numpy() # Convert to numpy array


        if fitIntercept:
            ones = np.ones((X.shape[0], 1))
            X = np.concatenate((ones, X), axis=1)

        n = X.shape[1] # Get number of covariates
        
        XY = np.column_stack((X, Y)) # Join X and Y

        # Calculate the SVD of XY
        _, _, VT_XY = np.linalg.svd(XY)
        V_XY = VT_XY.T
        Vxy = V_XY[0:n, n:] 
        Vyy = V_XY[n:, n:] 


        #Calculate the TLS estimator, and predictions
        betas = -np.divide(Vxy,Vyy) # The regression coefficients of TLS regression
        errorsXY = (-XY @ V_XY[:, n:] @ V_XY[:, n:].T) # The matrix of errors of X and Y
        errorX:np.ndarray = errorsXY[:, 0:n]
        errorY:np.ndarray = errorsXY[:, n:]
        predX = (X + errorX)
        predY = predX @ betas
        residuals = np.linalg.norm(errorsXY, axis=1)# Given by the frobenius Norm of the matrix of error of X and Y

        # print( f"residuals: \n {residuals}\n β_n:\n {betas}\n predictedY:\n {predY} \n predictedX:\n {predX} ")

        residuals = pd.DataFrame(residuals, index=samples, columns=['residualsTLS'])
        betas = pd.DataFrame(betas, index=features, columns=['betasTLS'])
        predY = pd.DataFrame(predY, index=samples, columns=['predYTLS'])
        predX = pd.DataFrame(predX, index=samples, columns=[f"{feature}_predXTLS"for feature in features])
        
        return residuals, betas, predY, predX
    

    def fit(self):
        """Fits the TLS regression model.

        Returns:
            self.residuals, self.betas, self.predY, self.predX
        """     
        self.residuals, self.betas, self.predY, self.predX = self.tlsRegression(self.X, self.Y)
        return self.residuals, self.betas, self.predY, self.predX   


    
    # def logLikelihoodTest(self):

    #     X = self.X.to_numpy() # Load Data
    #     Y = self.Y.to_numpy()
    #     M = self.M.to_numpy()

    #     # large Model

    #     if self.M2 is not None: # Use M2 if it exists as confounding factor
    #         M2 = self.M2.to_numpy()
    #         M = np.concatenate([M, M2], axis=1)

    #     largeModel = np.concatenate([M, X], axis=1) # Join all covariates

    #     # Small Model

    #     if self.M2 is not None: # Use M2 if it exists as confounding factor
    #         M2 = self.M2.to_numpy()
    #         M = np.concatenate([M, M2], axis=1)
    #     smallModel = M

    #     #TODO: Understand how to calculate the logLikelihoodTest when in TLS regression






class ResidualsLinearModel(GeneralLinearModel):

    def __init__(self, Y, X, M, M2=None, fit_intercept=True, copy_X=True, n_jobs=4, verbose=1, residualsType:str = "TLS"):

        assert residualsType in ["TLS", "malahanobis", None], "residualsType must be either TLS, None or malahanobis"
        
        if residualsType is not None:
            residualsCols = [col for col in X.columns if col[1] == residualsType]
            X = X.loc[:, residualsCols]
            X.columns = ['-'.join(col) for col in X.columns]
        super().__init__(Y, X, M, M2, fit_intercept, copy_X, n_jobs, verbose)


    def volcanoPlot(self, filepath:str, falseDiscoveryRate:float=0.01, pValHzLine:float = 0.001):
        """
        Volcano plot in order to find statisticall relevant relationships.
        """
        data = self.data.copy()
        data = data.loc[data['fdr'] < falseDiscoveryRate]
        yValues = -np.log10(data['pval'])
        xValues = data['beta']


        # Plot
        plt.scatter(
            xValues,
            yValues,
            c="k",
            s=5,
            alpha=0.5,
            rasterized=True,
        )

        # Labels
        plt.xlabel(r"$\beta$")
        plt.ylabel(r"$-\log_{10}(p-value)$")


        # Grid
        plt.axvline(0, c="k", lw=0.5, ls="--")
        plt.axhline(-np.log10(pValHzLine), c="k", lw=0.5, ls="--", label=f"p-value = {pValHzLine}")

        # Title
        plt.title("Volcano plot")
        plt.legend()
        self.volcanoPath = filepath
        plt.savefig(filepath, bbox_inches="tight")
        plt.close()

    def plotSignificantAssociation(self, proteomics:ProteinsMatrix, drugResponse:DrugResponseMatrix, filepath:str):

        results = self.data.copy()

        results = results.sort_values(by='pval', ascending=True)
        drug = results['y_id'].iloc[0]
        ppi  = results['x_id'].iloc[0]

        proteomics.plotPxPyDrug(drug, ppi, drugResponse,filepath)

    
class UnbiasedResidualsLinModel(MatrixData):
    """This class is a Linear Model that uses the residuals of the linear regression between X and Y as the new predictor for Drug response.
    But it does it by taking into consideration some confounding factors M, in our case growth props.
    """

    def __init__(self, ppis:Iterable[tuple[str,str]], proteomics:ProteinsMatrix, drugRes:DrugResponseMatrix, M:pd.DataFrame|pd.Series, fitIntercept=True, copy_X=True, standardisePx = True):
        """Will create linear model, that will be used for the two fits.

        Args:
            Y (_type_): Protein A expression
            X (_type_): Protein B expression
            M (_type_): Confounding Factors in Dataframe
            drugRes (DrugResponseMatrix): Drug Response Object
            fitIntercept (bool, optional): Does the Linear Model have an intercept parameter. Defaults to False.
            copy_X (bool, optional): _description_. Defaults to True.

        """
        self.ppis = ppis
        self.proteomics = proteomics
        self.M = M
        self.fitIntercept = fitIntercept
        self.copy_X = copy_X
        self.drugRes = drugRes
        self.drugRes.data =drugRes.data.T #Because in this object samples are columns
        self.standardisePx = standardisePx


    def checkCompleteObservations(self, X:pd.DataFrame, Y:pd.DataFrame):
        """ Removes samples that are not in common between X, Y, M and drug response matrices, which won't be useful

        Args:
            X (pd.DataFrame): Protein A expression
            Y (pd.DataFrame): Protein B expression

        Returns:
            _type_: (X, Y , self.M) with only common samples
        """ 
        X.dropna(axis=0, how='any', inplace=True)
        Y.dropna(axis=0, how='any', inplace=True)
        M :pd.DataFrame = self.M.dropna(axis=0, how='any')
        #Check if M is Categorical
        try:
            categoricalCols = M.select_dtypes(include=['object']).columns            
            M = pd.get_dummies(M, columns=categoricalCols)
            M.drop(columns=categoricalCols, inplace=True)
        except AttributeError as e:
            print(f"{e} M is a pd.Series")
            M = pd.get_dummies(M)
            M = M.iloc[:,1:]


        self.samples = set.intersection(set(X.index), set(Y.index), set(M.index), set(self.drugRes.data.index)) # Get common samples
        
        if len(self.samples) < 50:
            print(f"Not enough samples, {len(self.samples)}, in common between X, Y, M and drug response matrices \n Skipping")
            return None
        
        X = X.loc[self.samples]
        Y = Y.loc[self.samples]
        M = M.loc[self.samples]
        

        return X, Y, M


    

    def twoFits(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """ Performs two fits the first TLS Model with Py ~ M + Px and the second gLM with DrugResponse ~ M + residuals of the first fit

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: self.firstModelResiduals, self.firstModelBetas, self.secondModelResults
        """
        # Get the residuals of the first fit

        self.firstModelResiduals, self.firstModelBetas, M = self.fitPxPy() # This will give the coufounding factors M dummified

        # Fit the second model using the residuals of the first model as the new predictor for drug response
        Y = self.drugRes.data
        Y.fillna(Y.mean(), inplace=True)

        self.secondModel:GeneralLinearModel = GeneralLinearModel(self.drugRes.data, self.firstModelResiduals, M, fitIntercept=self.fitIntercept, copy_X=self.copy_X)

        self.secondModelResults = self.secondModel.fit_matrix()

        return self.firstModelResiduals, self.firstModelBetas, self.secondModelResults


    def fitPxPy(self):
        """ Fits the first model Py ~ M + Px, recursively for every PPI in self.ppis.
        It then returns the residuals of this recursion, for every PPI in self.ppis. Along with the coeficients of the regression.

        Returns:
           self.firstModelResiduals, self.firstModelBetas tuple[pd.DataFrame,pd.DataFrame]: Residuals of the first model and the coeficients of the regression
        """
        

        for XName, YName in self.ppis:
            X = self.proteomics.data.loc[:,XName]
            X = X.rename('proteinX')

            if self.standardisePx:
                X = (X - X.mean()) / X.std() # Standardise X

            Y = self.proteomics.data.loc[:,YName]

            XYMComplete = self.checkCompleteObservations(X, Y)
            if XYMComplete is None:
                print(f"X is {XName} and Y is {YName}, and there are no common samples between them, Drug Response and M \n Skipping")
                continue
            else:
                X, Y, M = XYMComplete
            
            X = pd.concat([X, M], axis=1)

            model:TLSRegression = TLSRegression(Y, X, fitIntercept=self.fitIntercept, copy_X=self.copy_X, standardise=False)
            residuals,betas,_,_ = model.fit()
            residuals.columns = [f"{YName};{XName}"] # Rename columns to match the PPI
            betas.columns = pd.MultiIndex.from_product([[YName], [XName]]) # Rename columns to match the PPI
            betas = betas.T
            residuals.index.name = "sample"
            betas.index.name = "sample"

            try:# If residuals already exists, merge them
                self.firstModelResiduals = pd.merge(self.firstModelResiduals, residuals, how='outer', on='sample')
            except AttributeError:
                self.firstModelResiduals = residuals

            try:# If betas already exists, concatenate them
                self.firstModelBetas = pd.concat([self.firstModelBetas, betas], axis=0)
            except AttributeError:
                self.firstModelBetas = betas



        assert self.firstModelResiduals is not None and self.firstModelBetas is not None, f"There are no Py ~ Px + M that has overlaping samples, so we can't build the first of the models"

        return self.firstModelResiduals, self.firstModelBetas, M
    

    def plotSignificantAssociation(self, proteomics:ProteinsMatrix, drugResponse:DrugResponseMatrix, filepath:str):

        results = self.secondModelResults.copy()

        results = results.sort_values(by='pval', ascending=True)
        drug = results['y_id'].iloc[0]
        ppi  = results['x_id'].iloc[0]
        if len(ppi.split(':')) > 0:
            ppi = ppi.split(':')[0] + ';' + ppi.split(':')[1]

        proteomics.plotPxPyDrug(drug, ppi, drugResponse, filepath)

def processPPIWrapper(self, ppi:tuple[str, str]) -> dict:
    """Wrapper for fitting the 2 linear models of Py ~ Px and Px ~ Py, so that it can be used in a multiprocessing pool

    Args:
        ppi (tuple[str, str]): Names of Py and Px
    Returns:
        dict: The results of the 2 linear models, one for Py ~ Px and the other for Px ~ Py
    """    
    results = [] # List of results for each drug
    for drugName in self.drugRes:
        
        YName = ppi[0]
        XName = ppi[1]
        _, _, res= self.getLinearModels(YName, XName, drugName)

        correctedPValues = multipletests(res[('info', 'logLikePValue')], method="fdr_bh")[1]
        res[('info', 'fdr')] = correctedPValues
        results.append(res)

        # invert Px and Py to understand if there are one way relationships
        YName = ppi[1]
        XName = ppi[0]
        _, _, res= self.getLinearModels(YName, XName, drugName)
        correctedPValues = multipletests(res[('info', 'logLikePValue')], method="fdr_bh")[1]
        res[('info', 'fdr')] = correctedPValues
        results.append(res)


    return results


class DRInteractionPxModel(MatrixData):
    """Linear Model Designed to find interactions between drug response and proteomics data, so the goal is to see what Drug Responses are impacted by a certain ppi (pY, pX)

    Args:
        MatrixData (_type_): _description_

    Returns:
        _type_: _description_
    """
    def __init__(self, ppis:Iterable[tuple[str,str]], proteomics:ProteinsMatrix, drugRes:DrugResponseMatrix, M:pd.DataFrame|pd.Series, fitIntercept=True, copyX=True, standardisePx = True, nJobs:int=4, filepath:str=None, data:pd.DataFrame=None, **readerKwargs):
        
        super().__init__(filepath, data, **readerKwargs)
        self.ppis = ppis
        self.proteomics = proteomics.data
        self.drugRes = drugRes.data
        self.M = M
        self.fitIntercept = fitIntercept
        self.copyX = copyX
        self.standardisePx = standardisePx
        self.nJobs = nJobs
        self.drugResLen = drugRes.data.shape[1]
        self.lenM =  M.shape[1]
    
    def modelRegressor(self):
        regressor = LinearRegression(
            fit_intercept=self.fitIntercept,
            copy_X=self.copyX,
            n_jobs=self.nJobs,
        )
        return regressor
    

    @staticmethod
    def loglike(y_true, y_pred):
        nobs = len(y_true)
        nobs2 = nobs / 2.0

        ssr = np.power(y_true - y_pred, 2).sum()

        llf = -nobs2 * np.log(2 * np.pi) - nobs2 * np.log(ssr / nobs) - nobs2

        return llf



    def getLinearModels(self, YName, XName, drugName) -> tuple[LinearRegression, LinearRegression, dict]:
        """Get the Linear Models (Larger and Smaller) for the given Protein X and Y Names
        It does this by subsetting the proteomics, drugRes, and M dataframes to only include samples common to all dataframes.
        And then builds the OLS models from statsmodels.api library.

        Args:
            YName (str): Protein Y name
            XName (str): Protein X name
            drugName (str): Drug name

        Returns:
            sm.OLS: larger Linear Model from statsmodels.api library
            sm.OLS: Smaller Linear Model from statsmodels.api library
            dict: Results of the linear model in dictionary format, with keys being effect size, p-value and other general info
        """
        
        Py = self.proteomics.loc[:,YName]
        Py = Py.dropna()
        Px = self.proteomics.loc[:,XName]
        Px = Px.dropna()
        M = self.M
        M = M.dropna(axis=0)
        drugRes = self.drugRes.loc[:,drugName]
        drugRes = drugRes.fillna(drugRes.mean())

        #get samples common to all dataframes
        samplesCommon = list(set.intersection(
            set(Py.index), set(Px.index), set(drugRes.index), set(M.index)
            ))# samples common to all dataframes
        
        #number of samples in common, n
        n = len(samplesCommon)


        #subset dataframes to common samples
        Py = Py.loc[samplesCommon]
        Px = Px.loc[samplesCommon]
        drugRes = drugRes.loc[samplesCommon]
        M = M.loc[samplesCommon]


        
        if self.standardisePx: # Zscore Px if standardisePx is True
            Px = (Px - Px.mean()) / Px.std()
            drugRes = (drugRes - drugRes.mean()) / drugRes.std()

        pxInteractionDR = drugRes.mul(Px, axis=0) # dR * Px


        #reordering of expressions to build the smaller and larger models
        # Small Model: Py ~ (Px + M) 
        # Large Model: Py ~ (Px + M) + (dr + Px:dR) 
        
        X = pd.concat([drugRes, pxInteractionDR], axis=1) 
        M = pd.concat([Px, M], axis=1)

        # Fit Confounding, small model
        lmSmall = self.modelRegressor().fit(M, Py)
        lmSmallLogLike = self.loglike(Py, lmSmall.predict(M))

        # Fit Confounding + features, Large model
        xLarge = pd.concat([M, X], axis=1)
        # Make sure all columns are strings
        xLarge.columns = xLarge.columns.astype(str)
        lmLarge = self.modelRegressor().fit(xLarge, Py)
        lmLargeLogLike = self.loglike(Py, lmLarge.predict(xLarge))

        # Log-ratio test
        lr = 2 * (lmLargeLogLike - lmSmallLogLike)
        LogLikeliRatioPVal = chi2(X.shape[1]).sf(lr)

        coefs = lmLarge.coef_
        columns = ['Px'] + self.M.columns.tolist() + ['drug'] + ['interaction']
        columns = [('effectSize', col) for col in columns]

        res = {col:coefs[index] for index,col in enumerate(columns)}
        res[('info', 'Py')] = YName
        res[('info', 'Px')] = XName
        res[('info', 'drug')] = drugName
        res[('info', 'n')] = n
        res[('info', 'logLikePValue')] = LogLikeliRatioPVal
        res[('info', 'llStatistic')] = lr
        res[('info', 'intercept')] = lmLarge.intercept_


        return lmLarge, lmSmall, res
    

    def fit(self)->pd.DataFrame:
        """Fit each Px and Py pair towards every drug in the drugRes dataframe.
            Calculate the Log likelihood p-value for the null hypothesis that the smaller model is correct, so the larger model does not add any covariate which is statistically significant.
            Or so to say the wilk's or likelihood ratio test.

        Returns:
            pd.DataFrame: The results of the fitting process, with the following columns: 
                Py, Px, drug, n, intercept, PxBeta, adherentBeta, semiAdherentBeta, suspensionBeta, unknownBeta, drugResBeta, interactionBeta, logLikePValue, llStatistic
        """        

        pararelList =  zip(repeat(self), self.ppis)
        with mp.Pool(CPUS) as process:
            pararelResults = process.starmap(processPPIWrapper, pararelList)
        results = list(chain.from_iterable(pararelResults))

        results = pd.DataFrame(results, columns = pd.MultiIndex.from_tuples(results[0].keys()))

        correctedPValues = multipletests(results['info']['logLikePValue'], method="fdr_bh")[1]
        results[('info', 'fdr')] = correctedPValues
        self.data = results

        return results
    
    def volcanoPlot(self, filepath:str, falseDiscoveryRate:float=0.01, pValHzLine:float = 0.001):
        """
        Volcano plot in order to find statisticall relevant relationships.
        """
        data = self.data.copy()
        data = data.loc[data['info']['fdr'] < falseDiscoveryRate]
        yValues = -np.log10(data['info']['logLikePValue'])
        xValues = data['effectSize']['interaction']


        plt.figure(figsize=(20, 20), dpi=600)
        # Plot
        plt.scatter(
            xValues,
            yValues,
            c="k",
            s=15,
            alpha=0.8,
            edgecolors="none",
            rasterized=True,
        )

        # Labels
        plt.xlabel(r"$\beta$")
        plt.ylabel(r"$-\log_{10}(p-value)$")

        # Grid
        plt.axvline(0, c="k", lw=0.5, ls="--")
        plt.axhline(-np.log10(pValHzLine), c="k", lw=0.5, ls="--", label=f"p-value = {pValHzLine}")

        # Title
        plt.title("Volcano plot")
        plt.legend()
        self.volcanoPath = filepath
        plt.savefig(filepath, bbox_inches="tight")
        plt.close()

        
    
    def scatterTheTopVolcano(self, filepathMold:str, proteomics:ProteinsMatrix, drugRes:DrugResponseMatrix, falseDiscoveryRate:float=0.10, topNumber:int=2):
        
        data = self.data.copy()
        data = data.loc[data['info']['fdr'] < falseDiscoveryRate]
        data = data.sort_values(by=[('info','logLikePValue'), ('effectSize','interaction')], ascending=[True, False])
        #Selecting top
        top = data.iloc[0:topNumber,:]
        #reseting index
        top = top.reset_index(drop=True)
        #iterate samples
        for index, row in top.iterrows():

            pValue = row['info']['logLikePValue']
            effectSize = row['effectSize']['interaction']
            drug = row['info']['drug']
            anotation = f'p-value: {pValue:.2e}\nβ: {effectSize:.2e} \ndrug: {drug} '
            anotation = {'text':anotation, 'xy':(0.1, 0.8), 'xycoords':'axes fraction', 'fontsize':10}
            filepath = filepathMold.split('.png')[0] + 'top'+ str(index) +'.png'
            ppi = row['info']['Py'] + ';' + row['info']['Px']
            proteomics.plotPxPyDrugContinous(drug, ppi, drugRes, filepath, **anotation)








    
    
        


    
