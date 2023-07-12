import pandas as pd
from resources import PATH, CPUS, read, DrugResponseMatrix, DRInteractionPxModel, ResiduesMatrix, ProteinsMatrix, PairwiseCorrMatrix
import statsmodels.api as sm
import statsmodels.formula.api as smf
from typing import Iterable
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.multitest import multipletests
import time as t

        

if __name__ == '__main__':

    #py ~ M + px + drugres + px:drugres

    drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    drugRes.data = drugRes.data.T
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    # using only corum ppis that we were able to recall, with high confidence
    vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')
    ppisOfInterest = set(vaeGLSPairwise.data.query("pValue < 0.00000001 & corum == 1").copy().index)
    ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}

    drugRes.data = drugRes.data.iloc[:, 0:1]
    ppisOfInterest = [('MRPL38', 'MRPL45')]

    M = pd.get_dummies(samplesheet['growth_properties'])
    M = M.rename(columns={'Semi-Adherent': 'SemiAdherent'})

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)


    dummy = DRInteractionPxModel(ppisOfInterest, ogProteomics, drugRes, M)
    dummy.fit()
    # dummy.write(PATH + '/internal/interactionModel/GLSPValueLessE-8VAEProteomics/regressor.pickle.gz')

