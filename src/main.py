# Imports
import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import *







if __name__ == '__main__':

    start  = t.time()
    baseModel :PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    proteomics: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    baseModel.data = baseModel.data.loc[baseModel.data['corum'] == 1]
    glsPairwise.data = glsPairwise.data.loc[glsPairwise.data['corum'] == 1]



    #PairwiseCorrMatrix.heatmap([baseModel, glsPairwise], ['pearsonR', 'beta'], [(0,1),(0,1)], 5, 'missingness', proteomics, 'heatmapMVperPPI5Bins.png')
    PairwiseCorrMatrix.heatmap([baseModel, glsPairwise], ['pearsonR', 'beta'], [(0,1),(0,1)], 10, 'missingness', proteomics, 'heatmapMVperPPI10Bins.png')
    

