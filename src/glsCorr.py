import numpy as np
import pandas as pd
from scipy.special import stdtr
import matplotlib.pyplot as plt
import time as t
from resources import *


if __name__ == '__main__':

    mofa2 = ProteinsMatrix(PATH + '/datasetsTese/proteomicsVAE.csv.gz', compression='gzip', index_col='Unnamed: 0')
    mofa2.write(PATH + '/datasetsTese/proteomicsVAE.pickle.gz')
    corum :ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
    
    pearsonPairwiseCorrs = mofa2.pearsonCorrelations('pearsonR') # get gls coefs and pearson r's for the inputed matrix
    glsPairwiseCorrs = mofa2.getGLSCorr()

    pearsonPairwiseCorrs.addGroundTruth(corum.ppis,'corum') #Add Ground truths and calculate AUC
    pearsonPairwiseCorrs.aucCalculator('corum', 'VAE-pearson AUC')

    glsPairwiseCorrs.addGroundTruth(corum.ppis,'corum') #Add Ground truths and calculate AUC
    glsPairwiseCorrs.aucCalculator('corum', 'VAE-gls AUC')

    pearsonPairwiseCorrs.write(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    glsPairwiseCorrs.write(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')

    drawRecallCurves([pearsonPairwiseCorrs, glsPairwiseCorrs],['blue', 'red'], '../images/VAEPearsonVsVAEGLSRecallCurve.png')


    
    

