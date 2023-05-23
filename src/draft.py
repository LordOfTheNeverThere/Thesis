# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *
import numpy as np


if __name__ == '__main__':


    corum:ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
    proteinMeanPairwiseCorr:PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelProteinMean.pickle.gz')
    print(proteinMeanPairwiseCorr.query('counts == 949 & corum == 1 & pValue == 0')) #18295 rows
    print(proteinMeanPairwiseCorr.query('counts == 949 & corum == 0 & pValue == 0')) #4408798 rows
    # proteinData :ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    # print(proteinData)
    # proteinData.data.dropna(axis=1, thresh=round(proteinData.data.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column
    # print(proteinData)
    # proteinData.data = proteinData.data.fillna(proteinData.data.mean())
    # print(proteinData)
    # mvPairwiseCorr : PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    # pValuePearson = mvPairwiseCorr.data['pValue']
    # glsPairwiseCorr : PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    # glsPairwiseCorr.aucCalculator('corum', 'gls Model', 'p-value', True)
    # pValueGLS = glsPairwiseCorr.data['p-value']
    # pairwiseCorr = proteinData.pearsonCorrelations('pearsonR')
    # pValuePearsonMean = pairwiseCorr.data['pValue']
    # pairwiseCorr.addGroundTruth(ppis=corum.ppis,externalDatasetName='corum')
    # pairwiseCorr.aucCalculator('corum', 'ProteinMean AUC', 'pValue', True)
    # glsPairwiseCorr.write(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    # pairwiseCorr.write(PATH + '/datasetsTese/baseModelProteinMean.pickle.gz')
    # drawRecallCurves([pairwiseCorr, glsPairwiseCorr],['blue', 'red'], '../images/pValuesOgMeanVsGLSRecallCurve.png')

    # plt.hist(pValuePearson, 30, alpha=0.5, label='p-values R')
    # plt.hist(pValuePearsonMean, 30, alpha=0.5, label='p-values R Mean')
    # plt.hist(pValueGLS, 30, alpha=0.5, label='p-values gls')

    # plt.xlabel('p-value')
    # plt.ylabel('Frequency of p-value')
    # plt.title('Histogram of p-values of different models')

    # plt.legend()

    # plt.savefig('../images/pValuesHist3Models.png')





    



    