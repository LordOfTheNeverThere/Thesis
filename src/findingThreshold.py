import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import time
from itertools import repeat
import multiprocessing as mp

from resources import *


RANDOMSTATE = None







# What Professor Pedro did not ask For :(

def getAUCvsThresholdPlot(pairwiseCorrData: pd.DataFrame) -> None:

    # find the best count threshold n
    best = True
    previousAUC = 0
    threshold = 2
    aucList = []
    thresholdList = []
    while threshold < 101:

        queriedPairCorrData = pairwiseCorrData.query('counts > @threshold')
        corrCumSum = np.cumsum(
            queriedPairCorrData['Corum']) / np.sum(queriedPairCorrData['Corum'])
        indexes = np.array(queriedPairCorrData.reset_index(
        ).index) / queriedPairCorrData.shape[0]
        currentAUC = auc(indexes, corrCumSum)
        # Update Lists for Curve plot
        aucList.append(currentAUC)
        thresholdList.append(threshold)

        if best and currentAUC < previousAUC:
            bestPairwuiseFrequency = threshold-1
            best = False
            print('The limit threshold of interactions were we start to lose information after increasing it is: ' +
                  str(threshold-1) + '\n with an auc of ' + str(previousAUC))
        # print(str(currentAUC) + '>=' + str(previousAUC) + '\n' + str(best))
        threshold += 1
        previousAUC = currentAUC

    # make chart with auc per n threshold

    _, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)
    ax.plot(
        thresholdList,
        aucList,
        c='green',
    )
    ax.plot([0, 100], [0.78, 0.78], "r--", label='AUC == 0.78', lw=0.7)
    ax.plot([bestPairwuiseFrequency, bestPairwuiseFrequency], [0, 1],
            "b--", label='max n==' + str(bestPairwuiseFrequency), lw=0.6)

    ax.legend(loc="lower right", frameon=False)
    ax.set_ybound(lower=0.6, upper=0.9)
    ax.set_ylabel("Area Under the Curve")
    ax.set_xlabel("Minimum interaction threshold")
    ax.set_title('AUC vs sample frequency of PPI threshold')
    ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

    plt.savefig("thresholdVSAUC.png",
                bbox_inches="tight")


# getAUCvsThresholdPlot(pairwiseCorrData)

    # What Professor Pedro Asked For :)

def variousRepeatsWrapper(iteration: int, sampleNum: int, proteinsData: ProteinsMatrix, glmCoefs: bool = False, pValueAUC:bool = False):
    corum:ppiDataset = read(PATH + '/external/ppiDataset/corum.pickle.gz')
    corum.name = 'corum'
    sampledProteins = ProteinsMatrix(None,  proteinsData.data.sample(n=sampleNum, axis=0, random_state=iteration * sampleNum))

    #get Pairwise Correlations
    if glmCoefs:
        pairwiseCorr = sampledProteins.getGLSCorr('Mean')
    else:
        pairwiseCorr = sampledProteins.pearsonCorrelations('coef', 'Mean', thresholdInteraction=3)
    #No longer need the proteinsData 
    del sampledProteins

    #Add all external Dataset
    pairwiseCorr.filepath = '.'
    PairwiseCorrMatrix.addGroundTruth(pairwiseCorr, corum.ppis, corum.name)
    # Calculate the AUC according to the proxyColumn, either p-value or coef
    if pValueAUC:
        AUC = pairwiseCorr.aucCalculator(corum.name, 'Mean', 'p-value', True)
    else:
        AUC = pairwiseCorr.aucCalculator(corum.name, 'Mean', 'coef', False)

    return AUC


def wrapperCheckPPIs(sampleNum: int, repeats: int, proteinsData: ProteinsMatrix, glmCoefs: bool = False, pValueAUC:bool = False):

    print(repeats, sampleNum)
    start = time.time()

    with mp.Pool(CPUS) as process:
        checkPPIGen = process.starmap(variousRepeatsWrapper, zip(range(0, repeats), repeat(sampleNum), repeat(proteinsData), repeat(glmCoefs), repeat(pValueAUC)))  # While Cycle
    result = list(checkPPIGen)
      
    print(time.time() - start)
    xLabel = str(sampleNum) +'| (n == ' + str(repeats) +')' 

    return xLabel, result


def randomSubSamplingAUC(proteinsData: ProteinsMatrix, subsampleSizes: list[int], repeats: list[int], glmCoefs:bool = False, pValueAUC:bool = False):


    checkPPIGen = map(wrapperCheckPPIs, subsampleSizes, repeats, repeat(
        proteinsData), repeat(glmCoefs), repeat(pValueAUC))  # First For Cycle

    allAUC = dict(checkPPIGen)

    if glmCoefs:
        globalPairwiseCorr: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')
    else:
        globalPairwiseCorr: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr75PV.pickle.gz')

    corumAUC = globalPairwiseCorr.aucs['p-value']['corum']
    _, ax = plt.subplots(figsize=(40, 8))
    ax.boxplot(allAUC.values(), labels=allAUC.keys())
    ax.set_ybound(lower=0.2, upper=1)
    ax.set_ylabel("AUC", fontsize=14)
    ax.set_xlabel("Sampling Number", fontsize=14)
    ax.axhline(y=corumAUC, color='red', linestyle='-', label='Corum AUC using p-values -> (Base Model))')
    ax.axhline(y=0.9*corumAUC, color='blue', linestyle=':', label="90% of Base Model's AUC")
    ax.axhline(y=0.8*corumAUC, color='blue',linestyle=':', label="80% of Base Model's AUC")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    plt.legend()
    if glmCoefs:
        plt.savefig("../images/glsAUCSubSampling_Mean_Proteomics.png",
                bbox_inches="tight")
    else:
        plt.savefig("../images/pearsonAUCSubSampling_Mean_Proteomics.png",
                bbox_inches="tight")

    plt.close()


if __name__ == '__main__':
    


    proteinsData: ProteinsMatrix = read(PATH + '/internal/proteomics/mean75PVProteomics.pickle.gz')
    corum:ppiDataset = read(PATH + '/external/ppiDataset/corum.pickle.gz')
    corum.name = 'corum'
    sampleNum = 185
    iterationNum = 1000
    aucList = []

    for iteration in range(0, iterationNum):
        start = time.time()
        sampledProteins = ProteinsMatrix(None,  proteinsData.data.sample(n=sampleNum, axis=0, random_state=RANDOMSTATE))
        pairwiseCorr = sampledProteins.getGLSCorr('Mean')
        PairwiseCorrMatrix.addGroundTruth(pairwiseCorr, corum.ppis, corum.name)
        aucList.append(pairwiseCorr.aucCalculator(corum.name, 'Mean', 'p-value', True))
        print(time.time() - start)

    aucData = pd.Series(aucList)

    globalPairwiseCorr: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')
    corumAUC = globalPairwiseCorr.aucs['p-value']['corum']
    _, ax = plt.subplots(figsize=(40, 8))

    ax.boxplot(aucData, labels=[f'{sampleNum}| (n == {iterationNum})'])
    ax.set_ylabel("AUC", fontsize=14)
    ax.set_xlabel("Sampling Number", fontsize=14)
    ax.axhline(y=corumAUC, color='red', linestyle='-', label='Corum AUC using p-values -> (Base Model))')
    ax.axhline(y=0.9*corumAUC, color='blue', linestyle=':', label="90% of Base Model's AUC")
    ax.axhline(y=0.8*corumAUC, color='blue',linestyle=':', label="80% of Base Model's AUC")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.show()



    # subsamplingList = list(range(5,950,5))
    # repeatsList= [round(900/repeat) + 5 if round(900/repeat) >= 4 and round(900/repeat) <= 100 else 100 if round(900/repeat)*2 > 100 else 5  for repeat in subsamplingList]

    # start = time.time()
    # randomSubSamplingAUC(proteinsData, subsamplingList, repeatsList, False, True)
    # end = time.time()
    # sumOfTime = end-start
    # print(sumOfTime)
    # start = time.time()
    # randomSubSamplingAUC(proteinsData, subsamplingList, repeatsList, True, True)
    # end = time.time()
    # sumOfTime = end-start
    # print(sumOfTime)

