
# Imports

import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix





if __name__ == '__main__':

    # get data to run lienar model with TLS residuals and plot significant associations with proteomics data
    drugRes = read(PATH + '/internal/drugResponse.pickle.gz')
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    vaeProteomics:ProteinsMatrix= read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz') 


    #Fit and test linear model for associations from VAE proteomics data's residuals
    residuals:ResiduesMatrix = read(PATH + '/internal/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')

    regressor = residuals.getLinearModel(drugRes, samplesheet)

    print(regressor.data['beta'].describe())
    regressor.plotSignificantAssociation(vaeProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSPValueLess0.001VAEProteomics.png')
    regressor.write(PATH + '/internal/residuals/GLSPValueLess0.001VAEProteomics/regressionMalahanobis.pickle.gz')

    #Fit and test linear model for associations from OG proteomics data's residuals

    pearsonResiduals = read(PATH + '/internal/residuals/pearsonPValueLess0.001OgProteomics/residuals.pickle.gz')

    pearsonrRegressor = pearsonResiduals.getLinearModel(drugRes, samplesheet)

    print(pearsonrRegressor.data['beta'].describe())
    pearsonrRegressor.plotSignificantAssociation(ogProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsPearsonPValueLess0.001OgProteomics.png')
    pearsonrRegressor.write(PATH + '/internal/residuals/pearsonPValueLess0.001OgProteomics/regressionMalahanobis.pickle.gz')
