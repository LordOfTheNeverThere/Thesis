
# Imports

import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix





if __name__ == '__main__':

    drugRes = read(PATH + '/datasetsTese/drugResponse.pickle.gz')
    samplesheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col=0)
    


    # regressionRes: ResidualsLinearModel = read(PATH + '/datasetsTese/regressionGlSGreater0.65MeanProteomics.pickle.gz')


    ogProteomics: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    # meanProteomics = read(PATH + '/datasetsTese/meanProteomics.pickle.gz') 
    vaeProteomics:ProteinsMatrix= read(PATH + '/datasetsTese/proteomicsVAE.pickle.gz') 


    # regressionRes.plotSignificantAssociation(meanProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSGreater0.65ProteinMean.png')^

    vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')

    # ppisOfInterest = set(vaeGLSPairwise.data.query("pValue < 0.001 & corum == 1").copy().index)
    # residuals = vaeProteomics.calculateResidues(ppisOfInterest)
    # residuals.write(PATH + '/datasetsTese/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')

    residuals:ResiduesMatrix = read(PATH + '/datasetsTese/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')

    regressor = residuals.getLinearModel(drugRes, samplesheet, 'malahanobis')

    print(regressor.data['beta'].describe())
    regressor.plotSignificantAssociation(vaeProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSPValueLess0.001VAEProteomics.png')
    regressor.write(PATH + '/datasetsTese/residuals/GLSPValueLess0.001VAEProteomics/regressionMalahanobis.pickle.gz')

### TODO


    baseModel :PairwiseCorrMatrix = read(PATH + '/internal/baseModelFiltered.pickle.gz')
    ppisOfInterest = set(baseModel.data.query("pValue < 0.001 & corum == 1").copy().index)
    print(len(ppisOfInterest))

    pearsonResiduals: ResiduesMatrix = ogProteomics.calculateResidues(ppisOfInterest)
    pearsonResiduals.write(PATH + '/internal/residuals/pearsonPValueLess0.001OgProteomics/residuals.pickle.gz')

    pearsonrRegressor = pearsonResiduals.getLinearModel(drugRes, samplesheet, 'malahanobis')

    print(pearsonrRegressor.data.iloc[:,-11:-1])
    pearsonrRegressor.plotSignificantAssociation(ogProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsPearsonPValueLess0.001OgProteomics.png')
    pearsonrRegressor.write(PATH + '/internal/residuals/pearsonPValueLess0.001OgProteomics/regressionMalahanobis.pickle.gz')
