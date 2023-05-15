# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *







if __name__ == '__main__':



    corum:ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
    proteinData :ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    proteinData.data.dropna(axis=1, thresh=round(proteinData.data.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column
    proteinData.data = proteinData.data.fillna(proteinData.data.mean())

    pairwiseCorr = proteinData.pearsonCorrelations('pearsonR')
    pairwiseCorr.addGroundTruth(ppis=corum.ppis,externalDatasetName='corum')
    pairwiseCorr.aucCalculator('corum', 'ProteinMean AUC')
    drawRecallCurves([pairwiseCorr, glsPairwise],['blue', 'red'], '../images/ogMeanVsGLSRecallCurve.png')







    # ogPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    # glsPairwise :PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    # ogProteomics :ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    # sampleSheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col='model_id')

    # # # # Plot the top 5 PPI that have a high pearson R and a low gls betas, and are not reported in corum, with points coloured by tissue type. 
    # # This will hopefully show a correlation due to tissue type and not due to the interaction of the members of the putative PPI. 
    # # Showing the protein baseline tecidual bias that pearson R cannot account for, a confounding factor.


    # # Merging/sorting and quering

    # mergedData:pd.DataFrame = ogPairwise.data.merge(right = glsPairwise.data[['glsCoefficient']], on ='PPI')
    # mergedData = mergedData.loc[mergedData['corum'] == 0]
    # mergedData = mergedData.loc[mergedData['biogrid'] == 0]
    # mergedData = mergedData.loc[mergedData['string'] == 0]
    # mergedData = mergedData.loc[mergedData['string400'] == 0]
    
    # mergedData = mergedData.query('globalCorrelation < 0.01 & globalCorrelation > 0').copy().sort_values(['glsCoefficient'], ascending=[False])
    
    # top5evidencePPI = mergedData.head(10)
    # top5evidencePPI.to_csv('../images/top10GlsbetasForHighR.csv')
    # indexes = list(top5evidencePPI.index)
    # setOfPPIs = [(proteins.split(';')[0], proteins.split(';')[1] ) for proteins in indexes]#Unpack PPIs of opposing inensities into tuples of proteins 
    # fig, ax = plt.subplots(10, 1, figsize=(15, 40))
    # for index, ppi in enumerate(setOfPPIs):

    #     proteinAB =  ogProteomics.data.loc[:,[ppi[0], ppi[1]]] #Get the protein expresion of the two proteins in the ppi
    #     proteinAB = proteinAB.dropna(axis=0)
    #     samples = list(proteinAB.index) #Get all samples where both proteins have an expression value
    #     tissues = [sampleSheet.loc[sample, 'tissue'] for sample in samples] #get all sample's tissues



    #     tissues = [TISSUEPALETTE[tissue] for tissue in tissues] #convert tissues to their respective colours


    #     ax[index].scatter(proteinAB[ppi[0]], proteinAB[ppi[1]], c = tissues)
    #     ax[index].set_xlabel(ppi[0])
    #     ax[index].set_ylabel(ppi[1])
    #     ax[index].set_title("Protein expression across samples, by Tissue, $r∈(0,0.001]$, $β_{GLS}↑$")
    #     ax[index].tick_params(labelsize=16)
    #     legend = [plt.Line2D([0], [0], marker='o', color=TISSUEPALETTE[key], label=key) for key in TISSUEPALETTE]
    #     ax[index].legend(handles = legend, fontsize=8, framealpha=0.2)

    # plt.savefig('../images/correlationsPearsonOnlyWithTissueTypev5.0.png')


