# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import utils


# %% Load Dataset

proteinsData = pd.read_csv('proteomicData.csv', index_col='modelId')

# %% Calculate Pearson Coorelation Value

# pearsonCorrRaw = proteinsData.corr(method='pearson')
# pearsonCorrRaw.to_csv('pearsonCorrRaw.csv', index=False)
pearsonCorrRaw = pd.read_csv('pearsonCorrRaw.csv')
# %% Get upwards triangular matrix data

pairwiseCorrRawSeries = utils.getPairwiseCorrData(data=pearsonCorrRaw)
pairwiseCorrRawSeries.sort_values(ascending=False, inplace=True) # Order it (descending Order)
pairwiseCorrRawSeries.to_csv('pairwiseCorrRawSeries.csv')
pairwiseCorrRawSeries =  pd.read_csv('pairwiseCorrRawSeries.csv', index_col='Unnamed: 0')

#%% Create columns for better integration with other external Databases by using the entrez ID of the subunits

pairwiseCorrRawData = utils.getGeneIDsCol(pairwiseCorrRawSeries)



# %% Load external Datasets

corumPPI = pd.read_json('../externalDatasets/corumPPI.json')
stringPPI = pd.read_table('../externalDatasets/stringPPI.txt', sep=' ')
geneIds = pd.read_table('../externalDatasets/geneIDsNCBI.tsv')

# %%
