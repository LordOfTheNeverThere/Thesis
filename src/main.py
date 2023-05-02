# Imports
import pandas as pd
import matplotlib.pyplot as plt
import random
from classes import ProteinsMatrix, PairwiseCorrMatrix, ppiDataset
from itertools import repeat
import utils

from env import PATH






if __name__ == '__main__':
    ogPairwise :PairwiseCorrMatrix = utils.read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = utils.read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    ogProteomics :ProteinsMatrix = utils.read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    sampleSheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col='model_id')

    # # # Plot the top 5 PPI that have a high pearson R and a low gls betas, and are not reported in corum, with points coloured by tissue type. 
    # This will hopefully show a correlation due to tissue type and not due to the interaction of the members of the putative PPI. 
    # Showing the protein baseline tecidual bias that pearson R cannot account for, a confounding factor.


    # Merging/sorting and quering

    mergedData:pd.DataFrame = ogPairwise.data[['globalCorrelation']].merge(right = glsPairwise.data[['glsCoefficient', 'corum']], on ='PPI')
    mergedData = mergedData.loc[mergedData['corum'] == 0]
    mergedData = mergedData.query('glsCoefficient < 0.2 & glsCoefficient > 0').copy().sort_values(['globalCorrelation'], ascending=[False])

    top5evidencePPI = mergedData.head(5)
    indexes = list(top5evidencePPI.index)
    setOfPPIs = [(proteins.split(';')[0], proteins.split(';')[1] ) for proteins in indexes]#Unpack PPIs of opposing inensities into tuples of proteins 
    fig, ax = plt.subplots(5, 1, figsize=(15, 40))
    for index, ppi in enumerate(setOfPPIs):

        proteinAB =  ogProteomics.data.loc[:,[ppi[0], ppi[1]]] #Get the protein expresion of the two proteins in the ppi
        proteinAB = proteinAB.dropna(axis=0)
        samples = list(proteinAB.index) #Get all samples where both proteins have an expression value
        tissues = [sampleSheet.loc[sample, 'tissue'] for sample in samples] #get all sample's tissues
        tissueColorDict = dict(zip(tissues,repeat('0')))

        usedColors = []
        for tissue in tissueColorDict:
            randColor = "#" + ''.join([random.choice('0123456789ABCDEF') for i in range(6)])
            while randColor in usedColors:
                randColor = "#" + ''.join([random.choice('0123456789ABCDEF') for i in range(6)])
            tissueColorDict[tissue] = randColor
            usedColors.append(randColor)


        tissues = [tissueColorDict[tissue] for tissue in tissues] #convert tissues to their respective colours


        ax[index].scatter(proteinAB[ppi[0]], proteinAB[ppi[1]], c = tissues)
        ax[index].set_xlabel(ppi[0])
        ax[index].set_ylabel(ppi[1])
        ax[index].set_title('Protein expression across samples, by Tissue')
        ax[index].tick_params(labelsize=16)
        legend = [plt.Line2D([0], [0], marker='o', color=tissueColorDict[key], label=key) for key in tissueColorDict]
        ax[index].legend(handles = legend, fontsize=8, framealpha=0.2)

    plt.savefig('../images/correlationsPearsonOnlyWithTissueType.png')