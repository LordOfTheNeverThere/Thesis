import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

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



# Proteomics PCA Scree and Cumulative plot


vaeProteomics = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')


pca = PCA(n_components=10).fit(vaeProteomics.data)

# Construct the plot of scree and the cumulative explained variance
# Get the explained variance and explained variance ratio
explained_variance = pca.explained_variance_
explained_variance_ratio = pca.explained_variance_ratio_

# Calculate the cumulative explained variance
cumulative_explained_variance = np.cumsum(explained_variance)

# Bar plot for explained variance ratios
ax1 = sns.barplot(x=np.arange(1, len(explained_variance_ratio) + 1), y=explained_variance_ratio, color='blue', alpha=0.8, edgecolor='k', linewidth=1, zorder=2)
ax1.set_xlabel('Principal Component')
ax1.set_ylabel('Explained Variance Ratio')
ax1.set_title('Scree Plot and Cumulative Explained Variance')
# Cumulative explained variance line plot    
ax2 = ax1.twinx()
ax2.plot(cumulative_explained_variance, marker='o', color='red')
ax2.set_xlabel('Cumulative Explained Variance')
ax2.grid(False)  
ax2.set_ylim(0)  
plt.tight_layout(pad=2)
plt.show()



# Good examples Scatter Plot of interactionModel without Drug in Small Model
dummy = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugLargeRegressor.pickle.gz')        

dummy.triangulate(0.175, 0.25, 27, 38, 13, 'goodExample.png', False)

# Good examples Scatter Plot of interactionModel with Drug in Small Model 
# first section
dummy = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugSmallRegressor.pickle.gz')        

triangulationResults = dummy.triangulate(0.099, 0.18, 50, 100, 16, 'goodExampleDrugSmall.png', False)

triangulationResults = dummy.triangulate(0.2, 0.3, 22, 50, 14, 'exampleDrugSmall0,2_0,3_22_50.png', True)