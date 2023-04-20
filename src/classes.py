from __future__ import annotations #  postpone evaluation of annotations
import pandas as pd
from itertools import combinations
import numpy as np
from sklearn.metrics import auc
from scipy.stats import pearsonr


class ppiDataset:

    def __init__(self, filepath, proteinLabels: list = [], **readerKwargs):

        self.data: pd.DataFrame = pd.read_csv(
            filepath, **readerKwargs)
        self.proteinLabels = proteinLabels
        self.ppis = set()

    def getPPIs(self, isCorum: bool = False) -> set:

        data = self.data.copy()

        if isCorum:

            def combinationsOfProteins(complx):

                if ';' in list(complx['subunits(Gene name)']):
                    complx['proteinTuple'] = list(combinations(
                        complx['subunits(Gene name)'].split(';'), 2))

                    return complx

            data = data.apply(
                lambda complx: combinationsOfProteins(complx), axis=1)
            ppiList = list(data.dropna()['proteinTuple'])
            ppiSet = {item for sublist in ppiList for item in sublist}

        else:
            data['proteinTuple'] = list(
                zip(data[self.proteinLabels[0]], data[self.proteinLabels[1]]))
            ppiSet = set(data['proteinTuple'])

        self.ppis = ppiSet

        return ppiSet


class ProteinsMatrix:

    def __init__(self, filepath: str, **readerKwargs):

        self.data: pd.DataFrame = pd.read_csv(
            filepath, **readerKwargs)


        

    def pearsonCorrelations(self, fileName: str, columnName: str, counting: bool = True, pValue: bool = True) -> pd.DataFrame:

        data = self.data.copy()
        # Correlation Matrix
        pearsonCorrMatrix = data.corr(method='pearson').round(decimals=4)
        print(pearsonCorrMatrix)

        pairwiseCorrData = pearsonCorrMatrix.where(np.triu(
            np.ones(pearsonCorrMatrix.shape), k=1).astype(bool)).stack().reset_index()
        print(pairwiseCorrData)
        pairwiseCorrData['PPI'] = pairwiseCorrData['level_0'] + \
            ";" + pairwiseCorrData['level_1']
        print(pairwiseCorrData)
        pairwiseCorrData.drop(columns=['level_0', 'level_1'], inplace=True)
        print(pairwiseCorrData)
        pairwiseCorrData = pairwiseCorrData.set_index('PPI')
        print(pairwiseCorrData)
        pairwiseCorrData.columns = [columnName]

        if counting:

            # Co-occorance Matrix
            # Get 1 and 0 depending on if there is a value or not in a specific spot
            coOccuranceMatrix = (data/data).fillna(0).astype(int)
            # Simple linear algebra to get the co-occurance values
            coOccuranceMatrix = coOccuranceMatrix.T.dot(coOccuranceMatrix)
            coOccuranceData = coOccuranceMatrix.where(np.triu(
                np.ones(coOccuranceMatrix.shape), k=1).astype(bool)).stack().reset_index()
            coOccuranceData['PPI'] = coOccuranceData['level_0'] + \
                ";" + coOccuranceData['level_1']
            coOccuranceData.drop(columns=['level_0', 'level_1'], inplace=True)
            coOccuranceData = coOccuranceData.set_index('PPI')
            coOccuranceData.columns = ['counts']
            pairwiseCorrData = pairwiseCorrData.merge(
                coOccuranceData, on='PPI', how='left')
            
        if pValue:

            def pearsonPValues(data:pd.DataFrame = None)-> pd.DataFrame|None:

                pValuesMatrix = data.corr(method=lambda x, y: pearsonr(x, y)[1]).round(decimals=3)
                pValuesMatrix = pValuesMatrix.where(np.triu(np.ones(pValuesMatrix.shape), k=1).astype(bool)).stack().reset_index()
                pairwisePValues = pValuesMatrix.drop(columns=['level_0', 'level_1'])
                pairwisePValues.columns = ['pValue']

                return pairwisePValues['pValue']
            
            pairwiseCorrData['pValue'] = pearsonPValues(pairwiseCorrData)
            print(pairwiseCorrData)
            

        pairwiseCorrData.sort_values(
            by=columnName, ascending=False, inplace=True)

        if fileName:
            pairwiseCorrData.to_csv(fileName + '.csv.gz', compression='gzip')

        self.data = pairwiseCorrData
        return pairwiseCorrData


class PairwiseCorrMatrix:

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, ** readerKwargs):


        if filepath:
            self.data: pd.DataFrame = pd.read_csv(filepath, **readerKwargs)
            
        elif data is not None:
            print('I am gooood!')
            self.data: data.copy()

        else:
            print('There should be either a filepath or data')
            return
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
