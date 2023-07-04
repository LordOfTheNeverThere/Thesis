
# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as t
from resources import ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix





if __name__ == '__main__':

    # get data to run lienar model with TLS residuals and plot significant associations with proteomics data
    drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    vaeProteomics:ProteinsMatrix= read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')


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


    #### Second Linear Model Py ~ M + Py and then D ~ M + residuals ####
    X = vaeProteomics.data.iloc[:,1].dropna().to_numpy()
    n = X.ndim
    Y = vaeProteomics.data.iloc[:,0].dropna().to_numpy()
    XY = np.vstack((X.T, Y)).T

    
    U_XY, S_XY, VT_XY = np.linalg.svd(XY)
    if n > 1:
        U_X, S_X, VT_X = np.linalg.svd(X)
        assert S_XY[-1] < S_X[-1]
    V_XY = VT_XY.T
    Vxy = V_XY[0:n, n:]
    Vyy = V_XY[n:, n:]

    betas = -Vxy / Vyy # The regression coefficients of TLS regression
    errorsXY = (-XY @ V_XY[:, n:] @ V_XY[:, n:].T) # The matrix of errors of X and Y
    errorX:np.ndarray = errorsXY[:, 0:n]
    errorY:np.ndarray = errorsXY[:, n:]
    predY = (X + errorX.T).T @ betas
    residuals = np.linalg.norm(errorsXY, axis=1)# Given by the frobenius Norm of the matrix of error of X and Y
    

    print( f"residuals: {residuals}\n β_n: {betas}\n predictedY: {predY} \n predictedX: {(X + errorX.T).T} ")