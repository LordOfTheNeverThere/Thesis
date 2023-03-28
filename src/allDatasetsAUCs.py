
# %%
import time as t
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from sklearn.metrics import auc
from classes import TreeNode

PATH = "../data"
# %%

pairwiseCorrData = pd.read_csv(
    PATH + '/datasetsTese/BaseModelPairwise.csv', index_col='PPI')

# # %%

# #Corum
# corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')
# #String
# stringPPI = pd.read_json(path_or_buf=PATH + 'stringPPI.json.gz', compression='gzip', index='String')
#biogrid
biogridPPI = pd.read_csv(PATH + '/externalDatasets/biogridPPIReduced.csv.gz', compression='gzip')

# import time as t




# listOfSets = list(map(set, biogridPPI.to_numpy()))

PPIs = TreeNode('root')
for PPI in biogridPPI.to_numpy():
    childrenValues = PPIs.getChildrenValue()

    if (PPI[0] in childrenValues) or (PPI[1] in childrenValues):

        if PPI[0] in childrenValues and PPI[1] not in childrenValues:
            proteinA = PPIs.getNodeFirstLayer(PPI[0])

            if PPI[1] not in proteinA.getChildrenValue():
                proteinA.addChild(TreeNode(PPI[1]))

        elif PPI[1] in childrenValues and PPI[0] not in childrenValues:

            proteinA = PPIs.getNodeFirstLayer(PPI[1])

            if PPI[0] not in proteinA.getChildrenValue():

                proteinA.addChild(TreeNode(PPI[0]))

    else:
        PPIs.addChild(TreeNode(PPI[0],{TreeNode(PPI[1])}))



    
print([[str(child1.value) + str(child.value)  for child1 in child.children] for child in PPIs.children])

# start = t.time()
# test = utils.addGroundTruth(listOfsets=listOfSets, data=pairwiseCorrData.head(2000), externalDatasetName='biogrid')
# end = t.time()
# print(end-start)
# print(test)            

