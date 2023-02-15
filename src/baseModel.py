# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt


#%% Load Dataset

proteinsData = pd.read_csv('proteomicData.csv', index_col='modelId')

# %% Calculate Pearson Coorelation Value

pearsonCoorRaw = protein
