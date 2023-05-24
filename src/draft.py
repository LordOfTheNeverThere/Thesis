# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *
import numpy as np


if __name__ == '__main__':

    limitsR = (-0.01,0.01)
    limitsGLS = (-0.5,0.5)
    corum = 1
    lung = True
    sortingColumn = 'glsCoefficient'
    ascending = False
    sortingSymbol = '↑' if ascending else '↓'


    baseModel :PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    ogProteomics :ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    sampleSheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col='model_id')

    # Merging/sorting and quering

    mergedData:pd.DataFrame = baseModel.data.merge(right = glsPairwise.data[['glsCoefficient']], on ='PPI')
    mergedData = mergedData.loc[mergedData['corum'] == corum]
    # mergedData = mergedData.loc[mergedData['biogrid'] == 0]
    # mergedData = mergedData.loc[mergedData['string'] == 0]
    # mergedData = mergedData.loc[mergedData['string400'] == 0]

    
    mergedData = mergedData.query('globalCorrelation < @limitsR[1] & globalCorrelation > @limitsR[0] & glsCoefficient > @limitsGLS[0] & glsCoefficient < @limitsGLS[1]').copy().sort_values([sortingColumn], ascending=[ascending])
    print(mergedData)
    
    top5evidencePPI = mergedData.head(5)
    indexes = list(top5evidencePPI.index)
    setOfPPIs = [(proteins.split(';')[0], proteins.split(';')[1] ) for proteins in indexes]#Unpack PPIs of opposing inensities into tuples of proteins 
    proteins =  [protein for sublist in setOfPPIs for protein in sublist]
    print(setOfPPIs, proteins)
    fig, ax = plt.subplots(10, 1, figsize=(20, 70))
    for index, protein in enumerate(proteins):
        
        proteinExpression =  ogProteomics.data.loc[:,protein] #Get the protein expresion of the two proteins in the ppi
        proteinExpression = proteinExpression.dropna(axis=0)

        samples = list(proteinExpression.index) #Get all samples where both proteins have an expression value
        mVSamples = [ogProteomics.data.loc[sample,:].isna().sum()  for sample in samples] #Count the #MV in each samples
        meanValsSamples = [ogProteomics.data.loc[sample,:].mean()  for sample in samples] #Count the #MV in each samples
        # replicateCountsSamples = [sampleSheet.loc[sample,'replicate_correlation']  for sample in samples] #Count the #MV in each samples
        # print(replicateCountsSamples)

        if not lung: # If we dont want a specific tisssue like lung

            samples = [sample for sample in samples if sampleSheet.loc[sample, 'tissue']!='Lung']
            proteinAB = proteinAB.loc[samples,:]
        
        tissues = [sampleSheet.loc[sample, 'tissue'] for sample in samples] #get all sample's tissues
        



        tissues = [TISSUEPALETTE[tissue] for tissue in tissues] #convert tissues to their respective colours


        ax[index].scatter(proteinExpression, mVSamples, c = tissues)
        ax[index].set_xlabel(protein)
        ax[index].set_ylabel('Mean of Sample')
        ax[index].set_title(f"Protein expression across samples, by Tissue, $r∈{limitsR}$, $β_{{GLS}}∈{limitsGLS}${sortingSymbol}, corum={corum}, lung={lung}")
        ax[index].tick_params(labelsize=16)
        legend = [plt.Line2D([0], [0], marker='o', color=TISSUEPALETTE[key], label=key) for key in TISSUEPALETTE]
        ax[index].legend(handles = legend, fontsize=8, framealpha=0.2)

    plt.savefig('../images/scatterPxvsMVTruePositives.png')






    



    