# Imports
import pandas as pd
from classes import ProteinsMatrix

from env import PATH

# Load Dataset

# proteinsData = pd.read_csv(
#     PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')




if __name__ == '__main__':

    proteinsData = ProteinsMatrix(PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col = 'modelID')
    proteinsData.pearsonCorrelations('baseModePairwiseWithpValues', 'globalCorrelation', counting=False)

