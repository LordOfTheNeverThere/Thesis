from itertools import combinations
from matrixData import MatrixData
import pandas as pd
import numpy as np



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


