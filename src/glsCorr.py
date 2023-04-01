import numpy as np
import pandas as pd
from scipy.special import stdtr
import matplotlib.pyplot as plt
import utils

PATH = "../data"

proteinData = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTrans.csv',index_col='modelID')

proteinData.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column

proteinDataMedian = proteinData.fillna(proteinData.median())
proteinDataMean = proteinData.fillna(proteinData.mean())


# proteinDataMean = proteinDataMean.drop(columns='modelID').to_numpy()
# proteinDataMedian = proteinDataMedian.drop(columns='modelID').to_numpy()


# Warp screen data and intercept based on covariance of screens

cholsigmainvMean = np.linalg.cholesky(np.linalg.inv(np.cov(proteinDataMean)))


warpedProteinsMean = proteinDataMean.T.values @ cholsigmainvMean
warpedIntereceptMean = cholsigmainvMean.T.sum(axis=0)


cholsigmainvMedian = np.linalg.cholesky(np.linalg.inv(np.cov(proteinDataMedian)))
warpedProteinsMedian = proteinDataMedian.T.values @ cholsigmainvMedian
warpedIntereceptMedian = cholsigmainvMedian.T.sum(axis=0)


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
GLS_p = 2 * stdtr(df, -np.abs(GLS_coef / GLS_se))
np.fill_diagonal(GLS_p, 1)

# # Save everything
# print(GLS_p)
# print(GLS_coef)

glsCorr = pd.DataFrame(GLS_coef)
pairwiseCorrData = glsCorr.where(np.triu(np.ones(glsCorr.shape), k=1).astype(bool)).stack().reset_index()
pairwiseCorrData.columns = ['index1', 'index2','glsBeta']
print(pairwiseCorrData)
pairwiseCorrData['glsBeta'].plot(kind='hist')

plt.savefig("testVariousAUC's.png", bbox_inches="tight")
np.save('GLS_p.npy', GLS_p)
np.save('GLS_sign.npy', np.sign(GLS_coef))
# screens.index.to_series().to_csv('genes.txt', index=False, header=False)
