# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *


if __name__ == '__main__':


    pearsonPairCorr:PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelPairwiseCorr.pickle.gz')
    pearsonPairCorrFiltered:PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')

    newFilteredData = pearsonPairCorrFiltered.data.merge(pearsonPairCorr.data['pValue'], how='left', on='PPI')
    newFilteredData['pValue_x'] = newFilteredData['pValue_y']
    newFilteredData.rename(columns={'pValue_x': 'pValue'}, inplace=True)
    newFilteredData.drop(columns=['pValue_y'], inplace=True)
    pearsonPairCorrFiltered.data = newFilteredData

    pearsonPairCorrFiltered.aucCalculator('corum', 'baseModel', 'pValue', ascending=True)
    pearsonPairCorr.aucCalculator('corum', 'baseModel', 'pValue', ascending=True)

    pearsonPairCorr.write(PATH + '/datasetsTese/baseModelPairwiseCorr.pickle.gz')
    pearsonPairCorrFiltered.write(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')




    