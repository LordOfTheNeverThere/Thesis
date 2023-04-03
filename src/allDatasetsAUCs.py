import time as t
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from sklearn.metrics import auc
from classes import TreeNode
import pickle

PATH = "../data"

pairwiseCorrData = pd.read_csv(PATH + '/datasetsTese/BaseModelPairwise.csv', index_col='PPI')
pairwiseCorrData = pairwiseCorrData.query('counts > 4') #Removing unmeaningful interactions


# #Corum
# corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')
# #String

# ppiStringFile = open('ppiStringTreeNode', 'rb')
# sringPPI = pickle.dump(ppiStringFile)
# ppiStringFile.close()

#biogrid
# ppiBiogridFile = open('ppiBiogridTreeNode', 'rb')
# biogridPPI = pickle.load(ppiBiogridFile)
# ppiBiogridFile.close()



# PPIs = TreeNode('root')
# for PPI in stringPPI.to_numpy():
#     childrenValues = PPIs.getChildrenValue()

#     if (PPI[0] in childrenValues) or (PPI[1] in childrenValues):

#         if PPI[0] in childrenValues and PPI[1] not in childrenValues:
#             proteinA = PPIs.getNodeFirstLayer(PPI[0])

#             if PPI[1] not in proteinA.getChildrenValue():
#                 proteinA.addChild(TreeNode(PPI[1]))

#         elif PPI[1] in childrenValues and PPI[0] not in childrenValues:

#             proteinA = PPIs.getNodeFirstLayer(PPI[1])

#             if PPI[0] not in proteinA.getChildrenValue():

#                 proteinA.addChild(TreeNode(PPI[0]))

#     else:
#         PPIs.addChild(TreeNode(PPI[0],{TreeNode(PPI[1])}))


# pairwiseCorrData = utils.addGroundTruthTreeNode(
#     ppiTree=biogridPPI, data=pairwiseCorrData, externalDatasetName="biogrid", filename='BaseModelPairwise')  # 150 sec


# pairwiseCorrData = utils.addGroundTruthTreeNode(
#     ppiTree=PPIs, data=pairwiseCorrData, externalDatasetName="string", filename='BaseModelPairwise')  # 150 sec





# Calculating Recall Curves

recallDict = dict()
indexes = np.array(pairwiseCorrData.reset_index().index) / pairwiseCorrData.shape[0]

for index, dbName in enumerate(['Corum', 'string', 'biogrid']):

    filterCorr = pairwiseCorrData[dbName]

    corrCumSum = np.cumsum(filterCorr) / np.sum(filterCorr)
    AUC = auc(indexes, corrCumSum)

    recallDict[dbName] = dict(x=list(indexes), y=list(
        corrCumSum), auc=AUC, color=sb.color_palette("tab20c").as_hex()[index])




_, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)

for ppiDB in recallDict:
    ax.plot(
    recallDict[ppiDB]['x'],
    recallDict[ppiDB]['y'],
    label="(AUC of " + str(ppiDB) + f"{recallDict[ppiDB]['auc']:.2f})",
        c=recallDict[ppiDB]['color'])



ax.plot([0, 1], [0, 1], "k--", lw=0.3)
ax.legend(loc="lower right", frameon=False)

ax.set_ylabel("Cumulative sum")
ax.set_xlabel("Ranked correlation")
ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

plt.savefig("combinedAllDatasetsRecallCurve.png",
            bbox_inches="tight")
plt.close("all")

_, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)

ax.barh(list(recallDict.keys()), [recallDict['Corum']['auc'], recallDict['string']['auc'], recallDict['biogrid']['auc']], height=0.2)

ax.set_xlabel("Area Under The Curve (AUC)")
ax.set_yticklabels(list(recallDict.keys()))
ax.set_xlim(left=0, right=1)
plt.savefig("combinedAllDatasetsBarplot.png",
            bbox_inches="tight")
plt.close("all")
