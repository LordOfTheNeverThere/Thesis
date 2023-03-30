import numpy as np
import pandas as pd
from scipy.special import stdtr
PATH = "../data"

proteinData = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTrans.csv',index_col='modelID')
proteinDataMedian = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTransMedian.csv.gz', compression='gzip')
proteinDataMean = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTransMean.csv.gz', compression='gzip')


proteinDataMean = proteinDataMean.drop(columns='modelID').to_numpy()
proteinDataMedian = proteinDataMedian.drop(columns='modelID').to_numpy()
proteinData.dropna(axis=1, inplace=True)

print(np.linalg.inv(np.cov(proteinData.T)))
# Warp screen data and intercept based on covariance of screens

cholsigmainvMean = np.linalg.cholesky(np.linalg.inv(np.cov(proteinDataMean.T)))
warpedProteinsMean = proteinDataMean.values @ cholsigmainvMean
warpedIntereceptMean = cholsigmainvMean.sum(axis=0)

cholsigmainvMedian = np.linalg.cholesky(np.linalg.inv(np.cov(proteinDataMedian.T)))
warpedProteinsMedian = proteinDataMedian.values @ cholsigmainvMedian
warpedIntereceptMedian = cholsigmainvMedian.sum(axis=0)


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
print(GLS_p)
print(GLS_coef)
np.save('GLS_p.npy', GLS_p)
np.save('GLS_sign.npy', np.sign(GLS_coef))
# screens.index.to_series().to_csv('genes.txt', index=False, header=False)
