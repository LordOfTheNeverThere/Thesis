# Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from resources import *
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import shapiro, pearsonr, chi2_contingency
from statsmodels.stats.diagnostic import het_white, het_breuschpagan


if __name__ == '__main__':


    proteomicsOG: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    vaeGLSPearson = read(PATH + '/internal/pairwiseCorrs/VAE/glsPearsonPairCorr.pickle.gz')
    setPPISOfInterest = list(vaeGLSPearson.data.sort_values(by='pValuePearson-GLS', ascending=False).query('corum == 1').head(5).index)

    featuresOfInterest = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0).columns[2:-1]


    for feature in ["age_at_sampling", "cancer_type", "cancer_type_detail", "tumour_grade", "sample_site", "mutational_burden" ,"cancer_type_ncit_id", "gender", "ethnicity", "growth_properties", "ploidy", "sampling_year", "supplier", "suppliers", "tissue", "tissue_status"]:
        print(feature)
        proteomicsOG.plotPxPySample(setPPISOfInterest, f'scatterPxPyAND{feature}.png', feature)

    # proteomics.write(PATH + '/internal/ProteinsMatrix/ogProteomics.pickle.gz')
    # vaeProteomics.write(PATH + '/internal/ProteinsMatrix/proteomicsVAE.pickle.gz')
    # meanProteomics.write(PATH + '/internal/ProteinsMatrix/meanProteomics.pickle.gz')
    # x = proteomics.iloc[:,10].dropna()
    # y = proteomics.iloc[:,15].dropna()
    # samplesCommon = x.index.intersection(y.index)
    # x = x.loc[samplesCommon]
    # X = pd.DataFrame({'x':x, 'intercept':np.ones(len(x))})
    # y = y.loc[samplesCommon]
    # pearsonr(x,y)

    # regressor = LinearRegression()
    # regressor.fit(x.values.reshape(-1,1), y.values.reshape(-1,1))
    # regressor.score(x.values.reshape(-1,1), y.values.reshape(-1,1))

    # residuals = y.values.reshape(-1,1) - regressor.predict(x.values.reshape(-1,1))
    # het_breuschpagan(residuals, X)
    # het_white(residuals, X)

    # plt.close()
    # plt.scatter(x, y)
    # plt.show()
    # plt.close()
    # plt.scatter(x, residuals)
    # plt.show()
