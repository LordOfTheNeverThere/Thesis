# Imports
from itertools import repeat
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
from resources import GeneDependency, DRPxPyInteractionPxModel, PyPxDrugInteractionModel, ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix, CPUS

def sampling(proteinsData: pd.DataFrame, sampleNum: int, iterationNum: int):
    corum = read(PATH + '/external/ppiDataset/corum.pickle.gz')
    corum.name = 'corum'
    sampledProteins = ProteinsMatrix(None,  proteinsData.sample(n=sampleNum, axis=0, random_state=iterationNum*sampleNum))
    pairwiseCorr = sampledProteins.getGLSCorr('Mean')
    PairwiseCorrMatrix.addGroundTruth(pairwiseCorr, corum.ppis, corum.name)
    auc = pairwiseCorr.aucCalculator(corum.name, 'Mean', 'p-value', True)
    return auc
if __name__ == '__main__':

    ogProteomics = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')

    #How was pv75MeanProteomics created?
    pv75MeanDf:pd.DataFrame = ogProteomics.data.copy()
    pv75MeanDf = pv75MeanDf.dropna(axis=1, thresh=round(0.75*pv75MeanDf.shape[0]))
    # How many nans per sample are there?
    sampleNans =  pv75MeanDf.isna().sum(axis=1)
    # sampleNans.hist(bins=100)
    # plt.show()
    # create set of 300 with most mvs
    missingSet = list(sampleNans.sort_values(ascending=False).head(350).index)
    # create set of 300 with least mvs
    presentSet = list(sampleNans.sort_values(ascending=True).head(350).index)
    sets = [missingSet, presentSet]

    iterationNum = 3
    sampleNums = [80, 250]
    proteinsData: ProteinsMatrix = read(PATH + '/internal/proteomics/mean75%PVProteomics.pickle.gz')
    globalPairwiseCorr: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')
    corumAUC = globalPairwiseCorr.aucs['p-value']['corum']

    
    aucData = {'auc':[], 'sampleNum':[], 'meanMVSet':[]}
    for sampleSet in sets:
        for sampleNum in sampleNums:
            
            pararelList =  zip(repeat(proteinsData.data.loc[sampleSet,:]), repeat(sampleNum), range(0, iterationNum))
            # use multiprocessing to speed up sampling
            
            with mp.Pool(CPUS) as process:
                aucList = process.starmap(sampling, pararelList)


            aucListLen = len(aucList)
            aucData['auc'] = aucData['auc'] + aucList
            aucData['sampleNum'] = aucData['sampleNum'] + list(repeat(sampleNum, aucListLen))
            aucData['meanMVSet'] = aucData['meanMVSet'] + list(repeat(round(sampleNans.loc[sampleSet].mean()), aucListLen))

    aucData = pd.DataFrame(aucData)
    grid = sns.FacetGrid(aucData, row="meanMVSet", col="sampleNum")

    grid.map(sns.boxplot, "auc")
    grid.set_axis_labels("auc", "")

    # add horizontal lines
    for ax in grid.axes.flat:
        ax.axvline(x=corumAUC, color='red', linestyle='-', label='Corum AUC using p-values -> (Base Model))')
        ax.axvline(x=0.9*corumAUC, color='blue', linestyle=':', label="90% of Base Model's AUC")
        ax.axvline(x=0.8*corumAUC, color='blue',linestyle=':', label="80% of Base Model's AUC")

    grid.add_legend()
    grid.set_titles(col_template="samplingNum={col_name}", row_template="MeanMVSamples{row_name}")
    plt.tight_layout(pad = 4)
    plt.subplots_adjust(hspace=0.5, wspace=0.2)
    plt.savefig("testingSubsampling.png", bbox_inches="tight")


