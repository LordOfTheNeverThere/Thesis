from matrixData import MatrixData
from resources import *
import math
import pandas as pd
import numpy as np  






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
        