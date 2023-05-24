from __future__ import annotations #  postpone evaluation of annotations
import pickle
import gzip
import pandas as pd
import numpy as np





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



class ResiduesMatrix(MatrixData):

    def __init__(self, filepath: str=None, data: pd.DataFrame=None, **readerKwargs):
        super().__init__(filepath, data, **readerKwargs)

