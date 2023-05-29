
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from resources import *



proteomics:ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
corum: ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
pearsonCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz') # R
vaeGLSCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz') # gls + vae
vaePearsonCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz') # R + vae


glsCorrs = proteomics.getGLSCorr(coefColumnName='beta')
glsCorrs.addGroundTruth(corum.ppis, 'corum') # gls

#Add all aucs

glsCorrs.aucCalculator('corum', 'gls','pValue', True)
glsCorrs.aucCalculator('corum', 'gls','beta', False)
print(glsCorrs)
glsCorrs.write(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')

pearsonCorrs.aucCalculator('corum', 'pearson','pValue', True)
pearsonCorrs.aucCalculator('corum', 'pearson','beta', False)
print(pearsonCorrs)
pearsonCorrs.write(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')

vaeGLSCorrs.aucCalculator('corum', 'VAE-GLS','pValue', True)
vaeGLSCorrs.aucCalculator('corum', 'VAE-GLS','beta', False)
print(vaeGLSCorrs)
vaeGLSCorrs.write(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')

vaePearsonCorrs.aucCalculator('corum', 'VAE-pearson','pValue', True)
vaePearsonCorrs.aucCalculator('corum', 'VAE-pearson','beta', False)
print(vaePearsonCorrs)
vaePearsonCorrs.write(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')





# tlsResidualsMatrix = proteomics.tlsResidues(glsPairwise)
# print(glsPairwise.data.loc[glsPairwise.data['p-value'] > 0])
# print(tlsResidualsMatrix.data.shape)
# print(tlsResidualsMatrix.data.isna().sum().sum())

