# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix

from env import PATH




if __name__ == '__main__':

    proteinsData = ProteinsMatrix(PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col = 'modelID')
    proteinsData.pearsonCorrelations('baseModePairwiseWithpValues', 'globalCorrelation', counting=False)

    pairwiseCorrPValues = PairwiseCorrMatrix(PATH + '/datasetsTese/baseModePairwiseWithpValues.csv.gz', compression = 'gzip')
    pairwiseCorr = PairwiseCorrMatrix(PATH + '/datasetsTese/BaseModelPairwise.csv')
    print(pairwiseCorrPValues.data)
    print(pairwiseCorr.data.merge(pairwiseCorrPValues.data, on='PPI', how='left'))
