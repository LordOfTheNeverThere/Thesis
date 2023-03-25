
# %%
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from sklearn.metrics import auc

PATH = "../data"
# %%

pairwiseCorrData = pd.read_csv(
    PATH + '/datasetsTese/BaseModelPairwise.csv', index_col='PPI')

# %%

#Corum
corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')
#String
stringPPI = pd.read_json(path_or_buf=PATH + 'stringPPI.json.gz', compression='gzip', index='String')
#biogrid
biogridPPI = pd.read_csv(PATH + '/externalDatasets/biogridPPIReduced.csv.gz', compression='gzip')

