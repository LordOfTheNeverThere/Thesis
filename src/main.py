# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as t
from statsmodels.stats.multitest import multipletests
from resources import UnbiasedResidualsLinModel, ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix




if __name__ == '__main__':



    # get gene dependecy data for effect size and its p-value

    effectSizes = pd.read_csv(PATH + '/internal/geneInteractionModel/CRISPRGeneEffectSize.csv', index_col=0)
    pValues = pd.read_csv(PATH + '/internal/geneInteractionModel/CRISPRGeneEffectSizePValue.csv', index_col=0)

    # get data to run lienar model with TLS residuals and plot significant associations with proteomics data
    # drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    # samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    # ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    # vaeProteomics:ProteinsMatrix= read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')
    # vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')



    # ppisOfInterest = set(vaeGLSPairwise.data.query("pValue < 0.001 & corum == 1").copy().index)
    # ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}

    
    # model1 = UnbiasedResidualsLinModel(ppisOfInterest, vaeProteomics, drugRes, samplesheet['growth_properties'], standardisePx=True)
    # fit = model1.twoFits()

    # model1:UnbiasedResidualsLinModel = read(PATH + '/internal/unbiasedResiduals/GLSPValueLess0.001VAEProteomics/regressor.pickle.gz')
    # model1.plotSignificantAssociation(ogProteomics, drugRes, 'unbiasedModelTest.png')

    #Fit and test linear model for associations from VAE proteomics data's residuals
    # residuals:ResiduesMatrix = read(PATH + '/internal/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')

    # regressor = residuals.getLinearModel(drugRes, samplesheet)

    # print(regressor.data.iloc[:, -10:-1].describe())
    # regressor.plotSignificantAssociation(vaeProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSPValueLess0.001VAEProteomics.png')
    # regressor.write(PATH + '/internal/residuals/GLSPValueLess0.001VAEProteomics/regressionMalahanobis.pickle.gz')

    #Fit and test linear model for associations from OG proteomics data's residuals

    # pearsonResiduals = read(PATH + '/internal/residuals/pearsonPValueLess0.001OgProteomics/residuals.pickle.gz')

    # pearsonrRegressor = pearsonResiduals.getLinearModel(drugRes, samplesheet)

    # print(pearsonrRegressor.data['beta'].describe())
    # pearsonrRegressor.plotSignificantAssociation(ogProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsPearsonPValueLess0.001OgProteomics.png')
    # pearsonrRegressor.write(PATH + '/internal/residuals/pearsonPValueLess0.001OgProteomics/regressionMalahanobis.pickle.gz')


