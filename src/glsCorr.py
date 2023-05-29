import numpy as np
import pandas as pd
from scipy.special import stdtr
import matplotlib.pyplot as plt
import time as t
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset


if __name__ == '__main__':

    proteomics: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')

    corum :ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
    
    glsPairwiseCorrs = proteomics.getGLSCorr()


    glsPairwiseCorrs.addGroundTruth(corum.ppis,'corum') #Add Ground truths and calculate AUC
    glsPairwiseCorrs.aucCalculator('corum', 'gls')

    pearsonPairwiseCorrs.write(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    glsPairwiseCorrs.write(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')

    drawRecallCurves([pearsonPairwiseCorrs, glsPairwiseCorrs],['blue', 'red'], '../images/VAEPearsonVsVAEGLSRecallCurve.png')


    
    

