# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import shapiro, pearsonr, chi2_contingency
from statsmodels.stats.diagnostic import het_white, het_breuschpagan


if __name__ == '__main__':

    # Create a sample matrix
    def isCovIdentity(data) -> None:

        covMatrix = abs(np.cov(data))

        # Create the expected identity matrix
        expected_matrix = np.eye(covMatrix.shape[0])

        # Perform the chi-square test
        chi2, p_value, _, _ = chi2_contingency(covMatrix)

        print("Chi-square statistic:", chi2)
        print("P-value:", p_value)


    proteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    vaeProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')
    meanProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/mean80%PVProteomics.pickle.gz')


    # proteomics.shapiroWilksTest()
    # vaeProteomics.shapiroWilksTest()
    # meanProteomics.shapiroWilksTest()
    isCovIdentity(proteomics.data)
    isCovIdentity(vaeProteomics.data)
    isCovIdentity(meanProteomics.data)

    vaeProteomics,_=ProteinsMatrix.whitening(vaeProteomics, vaeProteomics, saveIndexes=True)
    vaeProteomics = ProteinsMatrix(None, vaeProteomics)
    isCovIdentity(vaeProteomics.data)

    meanProteomics,_=ProteinsMatrix.whitening(meanProteomics, meanProteomics, saveIndexes=True)
    meanProteomics = ProteinsMatrix(None, meanProteomics)
    isCovIdentity(meanProteomics.data)
    # print(proteomics)


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
    



    