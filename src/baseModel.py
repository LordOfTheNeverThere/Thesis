# %% Imports
import pandas as pd
import seaborn as sb
# import numpy as np
# import matplotlib.pyplot as plt
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
pairwiseCorrRawData = pd.read_csv(
    '../data/datasetsTese/BaseModelPairwise.csv', index_col='Unnamed: 0')


# %% Load external Datasets

# corumPPI = pd.read_json('../data/externalDatasets/corumPPI.json')
# stringPPI = pd.read_table('../data/externalDatasets/stringPPI.txt', sep=' ')

# %% Dummy Test Pipeline function
_, _ = utils.getPairwiseCorrelation(fileName='Dummy', data=proteinsData.head(3), columnName='Test')

# %% Get get true or false on pairwise correlation
listOfSets = utils.getCorumListOfInteractions()
utils.addGroundTruth(listOfSets, pairwiseCorrRawData, 'Corum')


# %% Create Recall Curves

# %%
