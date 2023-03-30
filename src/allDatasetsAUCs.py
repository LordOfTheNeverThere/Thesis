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

pairwiseCorrData = pd.read_csv(
    PATH + '/datasetsTese/BaseModelPairwise.csv', index_col='PPI')


# #Corum
# corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')
# #String
stringPPI = pd.read_csv(PATH + '/externalDatasets/stringPPI900Selected.csv.gz', compression='gzip')
#biogrid
ppiBiogridFile = open('ppiBiogridTreeNode', 'rb')
biogridPPI = pickle.load(ppiBiogridFile)
ppiBiogridFile.close()

stringPPI.drop(columns=['String'],inplace=True)


PPIs = TreeNode('root')
for PPI in stringPPI.to_numpy():
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


pairwiseCorrData = utils.addGroundTruthTreeNode(
    ppiTree=biogridPPI, data=pairwiseCorrData, externalDatasetName="biogrid", filename='BaseModelPairwise')  # 150 sec

ppiStringFile = open('ppiStringTreeNode', 'ab')
pickle.dump(PPIs,ppiStringFile)
ppiStringFile.close()



pairwiseCorrData = utils.addGroundTruthTreeNode(
    ppiTree=PPIs, data=pairwiseCorrData, externalDatasetName="string", filename='BaseModelPairwise')  # 150 sec
