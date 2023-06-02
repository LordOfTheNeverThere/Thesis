# Imports
import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import *







if __name__ == '__main__':

    baseModel :PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    meanProteinsModel:PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelProteinMean.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    vaePearson: PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    vaeGLS: PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')
    
    barPlotCategories1 = {'$β_{GLS}$': {'ogMean \n coef': glsPairwise.aucs['beta'],'ogMean \n p-value': glsPairwise.aucs['pValue'], 'VAE \n coef': vaeGLS.aucs['beta'], 'VAE \n p-value': vaeGLS.aucs['pValue']},
                        'pearsonR': {'og \n coef': baseModel.aucs['pearsonR'],'og \n p-value': baseModel.aucs['pValue'], 'ogMean \n coef': meanProteinsModel.aucs['pearsonR'],'ogMean \n p-value': meanProteinsModel.aucs['pValue'], 'VAE \n coef': vaePearson.aucs['pearsonR'], 'VAE \n p-value': vaePearson.aucs['pValue']}
    }

    barPlotCategories2 = {'og \n pearsonR': {'coef': baseModel.aucs['pearsonR'], 'p-value': baseModel.aucs['pValue']},
                        'ogMean \n pearsonR': {'coef': meanProteinsModel.aucs['pearsonR'], 'p-value': meanProteinsModel.aucs['pValue']},
                        'ogMean \n $β_{GLS}$': {'coef': glsPairwise.aucs['beta'], 'p-value': glsPairwise.aucs['pValue']},
                        'VAE \n pearsonR': {'coef': vaePearson.aucs['pearsonR'], 'p-value': vaePearson.aucs['pValue']},
                        'VAE \n $β_{GLS}$': {'coef': vaeGLS.aucs['beta'], 'p-value': vaeGLS.aucs['pValue']},
    }
    


    plt.figure(figsize=(35, 35))
    
    pd.DataFrame(barPlotCategories2).plot(kind='bar')
    plt.tight_layout()
    plt.xticks(rotation = 0)
    plt.savefig('barplotAUCbyProxyMetric.png', bbox_inches="tight")
    pd.DataFrame(barPlotCategories1).plot(kind='bar')
    plt.tight_layout()
    plt.xticks(rotation = 0)
    plt.savefig('barplotAUCbyProteinsMatrix.png', bbox_inches="tight")


    # proteomics: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    # baseModel.data = baseModel.data.loc[baseModel.data['corum'] == 1]
    # glsPairwise.data = glsPairwise.data.loc[glsPairwise.data['corum'] == 1]



    # PairwiseCorrMatrix.heatmap([baseModel, glsPairwise], ['pearsonR', 'beta'], [(0,1),(0,1)], 5, 'missingness', proteomics, 'heatmapMVperPPI5Bins.png')
    # PairwiseCorrMatrix.heatmap([baseModel, glsPairwise], ['pearsonR', 'beta'], [(0,1),(0,1)], 10, 'missingness', proteomics, 'heatmapMVperPPI10Bins.png')
2