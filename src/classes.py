from __future__ import annotations #  postpone evaluation of annotations
import pandas as pd
from itertools import combinations
import numpy as np
from sklearn.metrics import auc
from scipy.stats import pearsonr
from env import PATH


class ppiDataset:

    def __init__(self, filepath, proteinLabels: list = [], **readerKwargs):

        self.data: pd.DataFrame = pd.read_csv(
            filepath, **readerKwargs)
        self.proteinLabels = proteinLabels
        self.ppis = set()

    def getPPIs(self, dataset:str) -> set:

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
            data = data.query("`Experimental System Type` == 'physical'")
            
            # Filter out Homedymers which are not object of our study

            data = data.query('`Official Symbol Interactor A` != `Official Symbol Interactor B`')

            data['proteinTuple'] = list(zip(data['Official Symbol Interactor A'], data['Official Symbol Interactor B']))
            ppiSet = set(data['proteinTuple'])

            self.ppis = ppiSet

        elif (dataset == allowedDatasets[3]):
            data['proteinTuple'] = list(zip(data['proteinA'], data['proteinB']))
            ppiSet = set(data['proteinTuple'])

            self.ppis = ppiSet



        return ppiSet


class ProteinsMatrix:

    def __init__(self, filepath: str, **readerKwargs):

        self.data: pd.DataFrame = pd.read_csv(
            filepath, **readerKwargs)


        

    def pearsonCorrelations(self, fileName: str, columnName: str, counting: bool = True, pValue: bool = True) -> PairwiseCorrMatrix:

        data = self.data.copy()
        # Get list with the names of every PPI
        proteinNames = data.columns.str.split(' ').str.get(0).to_numpy()
        ppiNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]
        # Correlation Matrix
        pearsonCorrMatrix = data.corr(method='pearson')

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

        if fileName:
            pairwiseCorrData.dropna().to_csv(PATH + '/datasetsTese/' + fileName + '.csv.gz', compression='gzip')

        return PairwiseCorrMatrix(None,pairwiseCorrData.dropna()) #There will be NAN correlations between proteins which do not appear simultaneously in at least two cell lines


class PairwiseCorrMatrix:

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, ** readerKwargs):


        self.data = data
        assert filepath or (
            data is not None), 'There should be either a filepath or data'

        if filepath:
            self.data: pd.DataFrame = pd.read_csv(filepath, **readerKwargs)
            
        elif data is not None:
            self.data: pd.DataFrame = data.copy()

        self.corrCumSum = None
        self.indexes = None
        self.auc = None
        self.label = None

    def addGroundTruth(self, ppis: set, externalDatasetName: str, filepath: str = None):
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

        if filepath:

            data.to_csv(filepath + '.csv.gz', compression='gzip')

        self.data = data

        return data
    
    def aucCalculator(self, yColumnName:str, label:str):

        pairwiseCorr = self.data 

        self.corrCumSum = np.cumsum(
            pairwiseCorr[yColumnName]) / np.sum(pairwiseCorr[yColumnName])
        
        self.indexes = np.array(pairwiseCorr.reset_index().index) / \
            pairwiseCorr.shape[0]
        self.auc = auc(self.indexes, self.corrCumSum)

        if not label:
            self.label = f"(AUC {self.auc:.2f})"
        
        self.label = label

    def query(self, query: str, inplace : bool = False) -> pd.DataFrame|None:

        if not inplace:
            return self.data.query(query)
        else:
            self.data = self.data.query(query)

    def compare(self, other: PairwiseCorrMatrix, querySelf: str, queryOther: str, key: str = 'PPI') -> pd.DataFrame:

        left: pd.DataFrame = self.query(querySelf)
        right: pd.DataFrame = other.query(queryOther)

        return left.merge(right, on=key)

        




        

# Deprecated
# class TreeNode:

#     def __init__(self, value, children: set = set()):

#         self.value = value
#         self.children = children

#     def addChild(self, childNode):

#         self.children.add(childNode)

#     def removeChild(self, childNode):

#         self.children = [
#             child for child in self.children if child is not childNode]

#     def transverse(self):

#         nodesToVisit = [self]

#         while len(nodesToVisit) > 0:

#             currentNode = nodesToVisit.pop()
#             print(currentNode.value)
#             # convert the set into a list so that the trasnverse works
#             nodesToVisit += list(currentNode.children)

#     def getNode(self, nodeValue):

#         if self.value == nodeValue:
#             return self

#         elif (len(self.children)):

#             for node in self.children:

#                 hasNode = node.getNode(nodeValue)
#                 if hasNode:
#                     return hasNode

#         return None

#     def getChildrenValue(self):
#         return {child.value for child in self.children}

#     def getNodeFirstLayer(self, nodeValue):

#         for child in self.children:

#             if child.value == nodeValue:

#                 return child
