# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix
import matplotlib.pyplot as plt

from env import PATH




if __name__ == '__main__':

    proteinsData = ProteinsMatrix(PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col = 'modelID')
    proteinsData.pearsonCorrelations('baseModePairwiseWithpValues', 'globalCorrelation')


    # pairwiseCorrPValues = PairwiseCorrMatrix(PATH + '/datasetsTese/baseModePairwiseWithpValues.csv.gz', compression = 'gzip')
    # ax = pairwiseCorrPValues.data.dropna()['pValue'].plot(kind='hist', bins=40, color='grey')
    # plt.savefig("../images/baseModelpValues.png",bbox_inches="tight")
    
