import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from sklearn.metrics import auc


PATH = "../data"

# Load Dataset

proteinsData = pd.read_csv(
    PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')



# Get Pairwise Correlation with n
pairwiseCorrData = utils.getPairwiseCorrelation(
    fileName='BaseModelPairwise', data=proteinsData, columnName='globalCorrelation')

# Load external Datasets
corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')


# Get get true or false on pairwise correlation
listOfSets = [set(subset.split(';')) for subset in corumPPI['subunits(Gene name)']]
pairwiseCorrData = utils.addGroundTruth(
    listOfSets, pairwiseCorrData, 'Corum', filename='BaseModelPairwise')



# find the best count threshold n
best = True
previousAUC = 0
threshold = 2
aucList=[]
thresholdList =[]
while threshold < 101:
    
    queriedPairCorrData = pairwiseCorrData.query('counts > @threshold')
    corrCumSum = np.cumsum(queriedPairCorrData['Corum']) / np.sum(queriedPairCorrData['Corum'])
    indexes = np.array(queriedPairCorrData.reset_index().index) / queriedPairCorrData.shape[0]
    currentAUC = auc(indexes, corrCumSum)
    #Update Lists for Curve plot
    aucList.append(currentAUC)
    thresholdList.append(threshold)
    
    
    if currentAUC < previousAUC:
        best = False
        print('The limit threshold of interactions were we start to lose information after increasing it is: ' + str(threshold-1) + '\n with an auc of ' + str(previousAUC))
    # print(str(currentAUC) + '>=' + str(previousAUC) + '\n' + str(best))
    threshold += 1
    previousAUC = currentAUC


# make chart with auc per n threshold

_, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)
ax.plot(
    thresholdList,
    aucList,
    label="BaseLine Model AUC with Various thresholds n",
    c='green',
)

ax.plot([0, 1], [0, 1], "k--", lw=0.3)
ax.legend(loc="lower right", frameon=False)

ax.set_ylabel("Area Under the Curve")
ax.set_xlabel("Minimum interaction threshold")
ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

plt.savefig("thresholdVSAUC.png",
            bbox_inches="tight")

