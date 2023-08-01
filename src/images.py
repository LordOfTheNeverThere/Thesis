import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt

from statsmodels.stats.multitest import multipletests
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset, PairwiseCorrMatrix



#How each image was created

# allPairwiseModels
#barplotAUCAllExternalDatasetsGLSPearson.png


aucs = pd.read_csv('aucData.csv', index_col=0)

#Build The Facet grid
# Matplotlib set main axis font size
plt.rcParams["axes.titlesize"] = 6

# Matplotlib set legend font size
plt.rcParams["legend.fontsize"] = 4

# Matplotlib set tick label font size
plt.rcParams["axes.labelsize"] = 4

# Matplotlib set tick label font size
plt.rcParams["xtick.labelsize"] = 4
plt.rcParams["ytick.labelsize"] = 4

grid = sns.FacetGrid(aucs, col="pairwiseDataset", row="metric")

grid.map(sns.barplot, "externalDataset", "auc", "gls", palette='viridis', errorbar=None)
#Add legend
grid.set_xticklabels(labels = ['Corum', 'Biogrid', 'String150', 'String400', 'String700', 'String900'] ,rotation=45)
grid.add_legend()
grid.set_titles(col_template="{col_name}", row_template="{row_name}")
grid.set_axis_labels("", "AUC")
plt.tight_layout(pad = 4)
plt.subplots_adjust(hspace=0.5, wspace=0.2)
plt.show()

