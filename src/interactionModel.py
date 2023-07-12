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
    ppisOfInterest = set(vaeGLSPairwise.data.query("corum == 1").sort_values(by='pValue', ascending=True).head(300).index)
    ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}



    M = pd.get_dummies(samplesheet['growth_properties'])
    M = M.rename(columns={'Semi-Adherent': 'SemiAdherent'})

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)


    # dummy = DRInteractionPxModel(ppisOfInterest, ogProteomics, drugRes, M)
    # fit = dummy.fit()
    # dummy.filepath = PATH + '/internal/interactionModel/GLSPValueVAEProteomicsHead300/regressor.pickle.gz'
    # dummy.write()

    dummy:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLSPValueVAEProteomicsHead300/regressor.pickle.gz')
    dummy.volcanoPlot('volcanoPlotDrInteractionPxModel.png') # 142393 points
    drugRes.data = drugRes.data.T
    dummy.scatterTheTop2Volcano('topVolcanoPlotScatter.png', ogProteomics, drugRes)

