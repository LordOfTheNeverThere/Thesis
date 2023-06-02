# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *







if __name__ == '__main__':

    
    baseModel :PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    proteomics: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')


    PairwiseCorrMatrix.heatmap([baseModel, glsPairwise], ['pearsonR', 'beta'], [(0,1),(0,1)], 10, 'missingness', proteomics)