# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils


# %% Load Dataset

proteinsData = pd.read_csv(
    '../data/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')


# %% Calculate Pearson Coorelation Value

# pearsonCorrRaw = proteinsData.corr(method='pearson')
# pearsonCorrRaw.to_csv('../data/datasetsTese/pearsonCorrRaw.csv', index=False)
# pearsonCorrRaw = pd.read_csv('../data/datasetsTese/pearsonCorrRaw.csv')
# %% Get upwards triangular matrix data

# pairwiseCorrRawSeries = utils.getPairwiseCorrData(data=pearsonCorrRaw)
# pairwiseCorrRawSeries.sort_values(ascending=False, inplace=True) # Order it (descending Order)
# pairwiseCorrRawSeries.to_csv('../data/datasetsTese/pairwiseCorrRawSeries.csv')
# pairwiseCorrRawSeries = pd.read_csv(
#     '../data/datasetsTese/pairwiseCorrRawSeries.csv', index_col='Unnamed: 0')

#%% Create columns for better integration with other external Databases by using the entrez ID of the subunits

# pairwiseCorrRawData = utils.getGeneIDsCol(pairwiseCorrRawSeries)
# pairwiseCorrRawData = pd.read_csv(
#     '../data/datasetsTese/BaseModelPairwise.csv', index_col='Unnamed: 0')


# %% Load external Datasets

corumPPI = pd.read_json('../data/externalDatasets/corumPPI.json')
# stringPPI = pd.read_table('../data/externalDatasets/stringPPI.txt', sep=' ')

# %% Pipeline function
# pairwiseCorrRawData = utils.getPairwiseCorrelation(
#     fileName='BaseModelPairwise', data=proteinsData, columnName='Correlation')
pairwiseCorrData = pd.read_csv(
    '../data/datasetsTese/BaseModelPairwise.csv', index_col='Unnamed: 0')

# %% Get get true or false on pairwise correlation
# listOfSets = [set(subset.split(';'))
#               for subset in corumPPI['subunits(Gene name)']]
# groundedPairwiseCorr = utils.addGroundTruth(
#     listOfSets, pairwiseCorrData, 'Corum', filename='BaseModelPairwise')


# %% Create Recall Curves
corrCumSum = np.cumsum(
    pairwiseCorrData['Corum']) / np.sum(pairwiseCorrData['Corum'])
indexes = np.array(pairwiseCorrData.reset_index().index) / \
    pairwiseCorrData.shape[0]
AUC = auc(indexes, corrCumSum)


_, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)
ax.plot(
    indexes,
    corrCumSum,
    label=f"(AUC {AUC:.2f})",
    c='blue',
)

ax.plot([0, 1], [0, 1], "k--", lw=0.3)
ax.legend(loc="lower right", frameon=False)

ax.set_ylabel("Cumulative sum")
ax.set_xlabel("Ranked correlation")
ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

# plt.savefig(f"{RPATH}/PPInteractions_roc_curves_overlap.pdf",
#             bbox_inches="tight")
plt.savefig("dummy.png",
            bbox_inches="tight")
plt.close("all")

# %%
