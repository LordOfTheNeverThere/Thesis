# Imports
import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import ResidualsLinearModel, read, PATH, ProteinsMatrix








if __name__ == '__main__':

    drugRes = read(PATH + '/datasetsTese/drugResponse.pickle.gz')
    


    regressionRes: ResidualsLinearModel = read(PATH + '/datasetsTese/regressionGlSGreater0.65MeanProteomics.pickle.gz')


    ogProteomics = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    meanProteomics = read(PATH + '/datasetsTese/meanProteomics.pickle.gz')  

    regressionRes.plotSignificantAssociation(meanProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSGreater0.65ProteinMean.png')