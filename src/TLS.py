import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as t
from statsmodels.stats.multitest import multipletests
from resources import GeneDependency, DRPxPyInteractionPxModel, ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix
from pathlib import Path

if __name__ == '__main__':

    drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    drugRes.data = drugRes.data.T
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    vaeProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz') #used for PCA computation
    ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz') #used for the interaction model class
           
    # Get ppis from String900 and Biogrid, to get ppis that are not ribossome enriched which is uninteresting
    vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')
    ppisOfInterest = set(vaeGLSPairwise.data.query("stringHighest == 1 & biogrid == 1 ").index)
    ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}

    pca, pcFactors = vaeProteomics.PCA(factorsName='PC', numPC=5)


    # Check if filepath exists
    folder = PATH + '/internal/interactionModel/String900andBiogrid/'
    assert Path(folder).is_dir(), 'Folder does not exist'

    dummy = DRPxPyInteractionPxModel(ppisOfInterest, ogProteomics, drugRes.data, pcFactors)
    #Testing writing object
    dummy.filepath = folder + 'interactionModelV.pickle.gz'
    dummy.write()
    start = t.time()
    fit = dummy.fit(numOfCores = 30)
    dummy.filepath = folder + 'interactionModelV.pickle.gz'
    dummy.write()
    print(f'fitting took {t.time() - start} seconds')

    





