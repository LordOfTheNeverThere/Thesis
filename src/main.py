
# Imports

import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix





if __name__ == '__main__':

    drugRes = read(PATH + '/datasetsTese/drugResponse.pickle.gz')
    samplesheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col=0)
    


    # regressionRes: ResidualsLinearModel = read(PATH + '/datasetsTese/regressionGlSGreater0.65MeanProteomics.pickle.gz')


    # ogProteomics = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    # meanProteomics = read(PATH + '/datasetsTese/meanProteomics.pickle.gz') 
    vaeProteomics:ProteinsMatrix= read(PATH + '/datasetsTese/proteomicsVAE.pickle.gz') 


    # regressionRes.plotSignificantAssociation(meanProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSGreater0.65ProteinMean.png')^

    vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')

    # ppisOfInterest = set(vaeGLSPairwise.data.query("pValue < 0.001 & corum == 1").copy().index)
    # residuals = vaeProteomics.calculateResidues(ppisOfInterest)
    # residuals.write(PATH + '/datasetsTese/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')

    residuals:ResiduesMatrix = read(PATH + '/datasetsTese/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')

    regressor = residuals.getLinearModel(drugRes, samplesheet, 'malahanobis')