import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from resources import PATH, read, DrugResponseMatrix, ResiduesMatrix, ProteinsMatrix, PairwiseCorrMatrix



if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)

    drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')

    residuals:ResiduesMatrix = read(PATH + '/internal/residuals/GLSPValueLess0.001VAEProteomics/residuals.pickle.gz')
    regressor = residuals.getLinearModel(drugRes, samplesheet)
    regressor.plotSignificantAssociation(ogProteomics, drugRes, 'pxpyTestingEffectSizeOfResidualsGLSPValueLess0.001VAEProteomics.png')
    regressor.write(PATH + '/internal/residuals/GLSPValueLess0.001VAEProteomics/regressor.pickle.gz')