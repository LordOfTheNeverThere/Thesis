import numpy as np
import pandas as pd
from scipy.special import stdtr
import matplotlib.pyplot as plt
import utils
from classes import PairwiseCorrMatrix, ProteinsMatrix

from env import PATH
def getGLSCorr(proteinData: ProteinsMatrix, pValues: bool = True) -> PairwiseCorrMatrix:

    
    proteinData = proteinData.data.copy()
    proteinData.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column

    proteinNames = proteinData.columns.str.split(' ').str.get(0).to_numpy()
    proteinNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]
    proteinDataMean = proteinData.fillna(proteinData.mean())

    cholsigmainvMean = np.linalg.cholesky(np.linalg.inv(np.cov(proteinDataMean)))
    warpedProteinsMean = proteinDataMean.T.values @ cholsigmainvMean
    warpedIntereceptMean = cholsigmainvMean.T.sum(axis=0)



    def linear_regression(warped_screens, warped_intercept):
        GLS_coef = np.empty((len(warped_screens), len(warped_screens)))
        GLS_se = np.empty((len(warped_screens), len(warped_screens)))
        ys = warped_screens.T
        for proteinIndex in range(len(warped_screens)):
            X = np.stack((warped_intercept, warped_screens[proteinIndex]), axis=1)
            coef, residues = np.linalg.lstsq(X, ys, rcond=None)[:2]
            df = warped_screens.shape[1] - 2
            GLS_coef[proteinIndex] = coef[1]
            GLS_se[proteinIndex] = \
                np.sqrt(np.linalg.pinv(X.T @ X)[1, 1] * residues / df)
        return GLS_coef, GLS_se


    GLS_coef, GLS_se = linear_regression(warpedProteinsMean, warpedIntereceptMean)
    df = warpedProteinsMean.shape[1] - 2
    #   Construct new PairwiseGLSCoefs Matrix

    glsCoefs = GLS_coef[np.triu_indices(GLS_coef.shape[0], k=1)]
    if pValues: #We might not want to add pValues so that we use less memory
        GLS_p = 2 * stdtr(df, -np.abs(GLS_coef / GLS_se))
        np.fill_diagonal(GLS_p, 1)
        glsPValues = GLS_p[np.triu_indices(GLS_p.shape[0], k=1)]
        pairwiseCorrData = pd.DataFrame(
            {'glsCoefficient': glsCoefs, 'p-value': glsPValues}, index=proteinNames)
    else:
        pairwiseCorrData = pd.DataFrame(
            {'glsCoefficient': glsCoefs}, index=proteinNames)
    
    pairwiseCorrData = PairwiseCorrMatrix(None, pairwiseCorrData)
    pairwiseCorrData.data.index.name = 'PPI'

    return pairwiseCorrData


if __name__ == '__main__':
    proteinData: ProteinsMatrix = utils.read(PATH + '/datasetsTese/ogProteomics.pickle.gz')

    getGLSCorr(proteinData)

