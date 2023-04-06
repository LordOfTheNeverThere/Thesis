# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from classes import ppiDataset

from env import PATH

# %% Load Dataset

proteinsData = pd.read_csv(
    PATH + '/datasetsTese/dummy.csv', index_col='modelID')

# %% Load external Datasets

corum = ppiDataset(filename=PATH + '/externalDatasets/corumPPI.csv.gz')
corum = corum.getPPIs(True)

pairwiseCorr = utils.getPairwiseCorrelation(
    fileName=None, data=proteinsData, columnName='mofaCorrelation')

pairwiseCorr = utils.addGroundTruth(corum, pairwiseCorr, 'Corum', None)

# %% Create Recall Curves
corrCumSum = np.cumsum(
    pairwiseCorr['Corum']) / np.sum(pairwiseCorr['Corum'])
indexes = np.array(pairwiseCorr.reset_index().index) / \
    pairwiseCorr.shape[0]
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
plt.savefig("mofaModelRecallCurve.png",
            bbox_inches="tight")
plt.close("all")
