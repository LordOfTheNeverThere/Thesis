# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *


if __name__ == '__main__':


    corum:ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
    proteinData :ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    proteinData.data.dropna(axis=1, thresh=round(proteinData.data.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column
    proteinData.data = proteinData.data.fillna(proteinData.data.mean())
    glsPairwiseCorr : PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    glsPairwiseCorr.aucCalculator('corum', 'gls Model', 'p-value', True)
    pairwiseCorr = proteinData.pearsonCorrelations('pearsonR')
    pairwiseCorr.addGroundTruth(ppis=corum.ppis,externalDatasetName='corum')
    pairwiseCorr.aucCalculator('corum', 'ProteinMean AUC')
    pairwiseCorr.write(PATH + '/datasetsTese/baseModelProteinMean.pickle.gz')
    drawRecallCurves([pairwiseCorr, glsPairwiseCorr],['blue', 'red'], '../images/ogMeanVsGLSRecallCurve.png')

    glsPairwiseCorr.aucCalculator()

    



    