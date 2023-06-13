import numpy as np
import pandas as pd
from scipy.special import stdtr
import matplotlib.pyplot as plt
import time as t
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset


if __name__ == '__main__':

    proteomics: ProteinsMatrix = read(PATH + '/internal/ogProteomics.pickle.gz')

    corum :ppiDataset = read(PATH + '/external/corum.pickle.gz')
    
    glsPairwiseCorrs = proteomics.getGLSCorr()


    glsPairwiseCorrs.addGroundTruth(corum.ppis,'corum') #Add Ground truths and calculate AUC
    glsPairwiseCorrs.aucCalculator('corum', 'gls')

    pearsonPairwiseCorrs.write(PATH + '/internal/VAEPearsonPairCorr.pickle.gz')
    glsPairwiseCorrs.write(PATH + '/internal/VAEGLSPairCorr.pickle.gz')

    drawRecallCurves([pearsonPairwiseCorrs, glsPairwiseCorrs],['blue', 'red'], '../images/VAEPearsonVsVAEGLSRecallCurve.png')


    
    

