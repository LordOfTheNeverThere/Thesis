# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix

from env import PATH




if __name__ == '__main__':

    proteinsData = ProteinsMatrix(PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col = 'modelID')
    proteinsData.pearsonCorrelations('baseModePairwiseWithpValues', 'globalCorrelation')
    # print(proteinsData.data)
    # pairwiseCorrPValues = PairwiseCorrMatrix(PATH + '/datasetsTese/baseModePairwiseWithpValues.csv.gz', compression = 'gzip')
    # print(pairwiseCorrPValues.data)
    
