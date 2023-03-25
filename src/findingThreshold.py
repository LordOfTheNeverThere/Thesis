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
while best:

    pairwiseCorrData.query('counts > @threshold')
    corrCumSum = np.cumsum( pairwiseCorrData['Corum']) / np.sum(pairwiseCorrData['Corum'])
    indexes = np.array(pairwiseCorrData.reset_index().index) / pairwiseCorrData.shape[0]
    currentAUC = auc(indexes, corrCumSum)
    
    if currentAUC < previousAUC:
        best = False
        print('The limit threshold of interactions were we start to lose information after increasing it is: ' + str(threshold-1) + '\n with an auc of ' + str(previousAUC))

    threshold += 1
    previousAUC = currentAUC
    

