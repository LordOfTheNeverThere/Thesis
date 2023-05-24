# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *







if __name__ == '__main__':
    limitsR = (0.9,1)
    limitsGLS = (-0.001,0.001)
    corum = 0
    lung = False
    sortingColumn = 'glsCoefficient'
    sortingType = False


    vaePairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    vaeGlsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')
    ogProteomics :ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    sampleSheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col='model_id')

    # # # Plot the top 5 PPI that have a high pearson R and a low gls betas, and are not reported in corum, with points coloured by tissue type. 
    # This will hopefully show a correlation due to tissue type and not due to the interaction of the members of the putative PPI. 
    # Showing the protein baseline tecidual bias that pearson R cannot account for, a confounding factor.


    # Merging/sorting and quering

    mergedData:pd.DataFrame = vaePairwise.data.merge(right = vaeGlsPairwise.data[['glsCoefficient']], on ='PPI')
    mergedData = mergedData.loc[mergedData['corum'] == corum]
    # mergedData = mergedData.loc[mergedData['biogrid'] == 0]
    # mergedData = mergedData.loc[mergedData['string'] == 0]
    # mergedData = mergedData.loc[mergedData['string400'] == 0]
    
    mergedData = mergedData.query('pearsonR < @limitsR[1] & pearsonR > @limitsR[0] & glsCoefficient > @limitsGLS[0] & glsCoefficient < @limitsGLS[1]').copy().sort_values([sortingColumn], ascending=[sortingType])
    print(mergedData)
    
    top5evidencePPI = mergedData.head(10)
    indexes = list(top5evidencePPI.index)
    setOfPPIs = [(proteins.split(';')[0], proteins.split(';')[1] ) for proteins in indexes]#Unpack PPIs of opposing inensities into tuples of proteins 
    fig, ax = plt.subplots(10, 1, figsize=(20, 70))
    for index, ppi in enumerate(setOfPPIs):

        proteinAB =  ogProteomics.data.loc[:,[ppi[0], ppi[1]]] #Get the protein expresion of the two proteins in the ppi
        proteinAB = proteinAB.dropna(axis=0)

        samples = list(proteinAB.index) #Get all samples where both proteins have an expression value

        if not lung: # If we dont want a specific tisssue like lung

            samples = [sample for sample in samples if sampleSheet.loc[sample, 'tissue']!='Lung']
            proteinAB = proteinAB.loc[samples,:]
        
        tissues = [sampleSheet.loc[sample, 'tissue'] for sample in samples] #get all sample's tissues
        



        tissues = [TISSUEPALETTE[tissue] for tissue in tissues] #convert tissues to their respective colours


        ax[index].scatter(proteinAB[ppi[0]], proteinAB[ppi[1]], c = tissues)
        ax[index].set_xlabel(ppi[0])
        ax[index].set_ylabel(ppi[1])
        ax[index].set_title(f"Protein expression across samples, by Tissue, $r∈{limitsR}$, $β_{{GLS}}∈{limitsGLS}$, corum={corum}, lung={lung}")
        ax[index].tick_params(labelsize=16)
        legend = [plt.Line2D([0], [0], marker='o', color=TISSUEPALETTE[key], label=key) for key in TISSUEPALETTE]
        ax[index].legend(handles = legend, fontsize=8, framealpha=0.2)

    plt.savefig('../images/correlationsPearsonOnlyWithTissueTypev8.1.png')


