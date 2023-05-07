#%% Imports 
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from resources import *
# %% Load file and Transform (Only Once)

# proteinsData = pd.read_csv('../datasetsTese/proteomics.csv', index_col='GeneSymbol')
# display(proteinsData)
# proteinsData = pd.DataFrame.transpose(proteinsData)
# proteinsData.index.names = ['modelId']
# display(proteinsData)
# proteinsData.to_csv('proteomicData.csv')

# drugResponseData = pd.read_csv('../datasetsTese/drugresponse.csv', index_col= 'drug')
# display(drugResponseData)
# drugResponseData = pd.DataFrame.transpose(drugResponseData)
# drugResponseData.index.names = ['modelId']
# display(drugResponseData)
# drugResponseData.to_csv('drugResponse.csv')

# %% Load Transformed File

proteinsData = pd.read_csv('proteomicData.csv', index_col='modelId')
drugResponseData = pd.read_csv('drugResponse.csv', index_col='modelId')

# %% Missing Value

fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle(
    'Number of Proteins with a certain Absence Throughout Samples')

numSamples = len(proteinsData.index.values) # Get num of samples so we can calculate the percentage of mv's
proteinsNAN = proteinsData.isna()
proteinsNAN = round((proteinsNAN.sum()/numSamples) * 100, 2)


axes = sb.histplot(data= proteinsNAN, bins=20)
axes.set(xlabel = 'Percentage of missing values', ylabel = 'Num Of Proteins')

fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle(
    'Number of Drugs with a certain Absence Throughout Samples')

# Get num of samples so we can calculate the percentage of mv's
numSamples = len(drugResponseData.index.values)
drugResponseNAN = drugResponseData.isna()
drugResponsesNAN = round((drugResponseNAN.sum()/numSamples) * 100, 2)


axes = sb.histplot(data=drugResponsesNAN, bins=20)
axes.set(xlabel='Percentage of missing values', ylabel='Num Of Drugs')

## The lack of presence of a Drug Response in the sample space in the case of a specific drug is not important, so no missing values treatement will be done.
## Outlier treatments in any of the datasets is not at all relevant, since a different value might due to somo important biomolecular cues...

# %% Segmentation by Mvs



truncatedProteinsNAN = proteinsNAN.where(proteinsNAN<=35) # Only count proteins with more than 65% incidence in samples
truncatedProteinsNAN = truncatedProteinsNAN.dropna()
relevantProteins = list(truncatedProteinsNAN.index)

relevantProteinData = proteinsData[relevantProteins]
# %% Distros


fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle('BoxPlots of Proteins with Incidence Greater than 65%')

plt.xticks(rotation = 90)

sb.histplot(data=relevantProteinData.head(10), ax=ax1)

# %% Data Description
stats = proteinsData.describe()
statsRelevant = relevantProteinData.describe()
statsDrugResponse = drugResponseData.describe()

# %% IQR Data Distribution
thirdQuartile = stats.iloc[-2]
firstQuartile = stats.iloc[-4]
iqrProteins = thirdQuartile - firstQuartile

thirdQuartile = statsRelevant.iloc[-2]
firstQuartile = statsRelevant.iloc[-4]
iqrRelevantProteins = thirdQuartile - firstQuartile

thirdQuartile = statsDrugResponse.iloc[-2]
firstQuartile = statsDrugResponse.iloc[-4]
iqrDrugResponse = thirdQuartile - firstQuartile

 #Join both in a dataframe
iqrConcatProteins = pd.concat([iqrProteins, iqrRelevantProteins], axis=1)
iqrConcatProteins.columns = ['IQR All Proteins',
                             'IQR Missing Values Truncated Proteins']

# %% Histogram of IQR's

fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle(
    'Inter-Quantile Range Distribution')


axes = sb.histplot(data=iqrConcatProteins, ax=ax1, bins=100)
axes.set(xlabel='Inter-Quantile Range', ylabel='Num Of Proteins')

fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
fig.suptitle(
    'Inter-Quantile Range Distribution')


axes = sb.histplot(data=iqrDrugResponse, ax=ax1, bins=100)
axes.set(xlabel='Inter-Quantile Range', ylabel='Num Of Drugs')


# %% Proteins that are prevalent and vary alot
truncatedIQRProteins = iqrRelevantProteins.where(iqrRelevantProteins > 2.5)
truncatedIQRProteins = truncatedIQRProteins.dropna()
truncatedIQRProteins.sort_index(inplace = True)
print(round(truncatedIQRProteins,2).to_latex())

truncatedIQRDrugs = iqrDrugResponse.where(iqrDrugResponse > 3)
truncatedIQRDrugs = truncatedIQRDrugs.dropna()
truncatedIQRDrugs.sort_values(inplace=True)
print(round(truncatedIQRDrugs, 2).to_latex())

# %%
