from matrixData import MatrixData
from sklearn.metrics import auc
import pandas as pd
import numpy as np
from resources import *


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
    
    def aucCalculator(self, yColumnName:str, label:str, proxyColumn:str, ascending:bool ):
        """Adds the value of AUC of the Recall curve using a specified external PPI dataset with yColumnName

        Args:
            yColumnName (str): Name of the df column where there is the truth value of the existence or not of the PPI in the reported external PPI dataset
            label (str): Text which will show up as label next to the value of the AUC, e.g 'Baseline Auc == 0.9' 
            proxyColumn (str): Name of the column of the statistical meausure to quantify the probability of PPI existence
            ascending(bool): should the proxyColumn be ordered ascending or not for best AUC calculation
        """
        pairwiseCorr = self.data 

        
        pairwiseCorr.sort_values(by=proxyColumn, ascending=ascending, inplace=True) # We sort rows by the smallest to greatest pValues

        self.corrCumSum = np.cumsum(
            pairwiseCorr[yColumnName]) / np.sum(pairwiseCorr[yColumnName])
        
        self.indexes = np.array(pairwiseCorr.reset_index().index) / \
            pairwiseCorr.shape[0]
        self.auc = auc(self.indexes, self.corrCumSum)

        if not label:
            self.label = f"(AUC {self.auc:.2f})"
        
        self.label = label + f" (AUC {self.auc:.2f})"

  