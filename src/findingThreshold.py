import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from sklearn.metrics import auc
from classes import ppiDataset
import time
from itertools import repeat
from env import PATH
import multiprocessing as mp


RANDOMSTATE = None
CPUS = 30
assert CPUS < mp.cpu_count() - 1


proteinsData = pd.read_csv(
    PATH+'/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')

globalPairwiseCorr = pd.read_csv(
    PATH+'/datasetsTese/BaseModelPairwise.csv', index_col='PPI')
globalPairwiseCorr = globalPairwiseCorr['Corum']




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

def variousRepeatsWrapper(iteration: int, sampleNum: int, proteinsData: pd.DataFrame, corum: pd.DataFrame):
    

    proteinsDatSample = proteinsData.sample(n=sampleNum, axis=0, random_state=iteration * sampleNum)
    pairwiseCorr = utils.getPairwiseCorrelation(proteinsDatSample, None, str(sampleNum) + ' samples', False)
    del proteinsDatSample
    pairwiseCorr = pairwiseCorr.merge(globalPairwiseCorr, on='PPI', how='left')
    # pairwiseCorr = utils.addGroundTruth(corum, pairwiseCorr.head(2000), 'Corum', None) DEPRECATED
    corrCumSum = np.cumsum(
        pairwiseCorr['Corum']) / np.sum(pairwiseCorr['Corum'])
    indexes = np.array(pairwiseCorr.reset_index().index) / \
        pairwiseCorr.shape[0]
    del pairwiseCorr
    AUC = auc(indexes, corrCumSum)
    return AUC


def wrapperCheckPPIs(sampleNum: int, repeats: int, proteinsData: pd.DataFrame, corum: pd.DataFrame):

    print(repeats, sampleNum)
    start = time.time()

    with mp.Pool(CPUS) as process:
        checkPPIGen = process.starmap(variousRepeatsWrapper, zip(range(0, repeats), repeat(
            sampleNum), repeat(proteinsData), repeat(corum)))  # While Cycle
    result = list(checkPPIGen)
    if sampleNum == 450:
        print('50% Done')
        
    print(time.time() - start)

    return sampleNum, result


def randomSubSamplingAUC(proteinsData: pd.DataFrame, subsampleSizes: list[int], repeats: list[int]):

    corum = ppiDataset(filename=PATH + '/externalDatasets/corumPPI.csv.gz')
    corum = corum.getPPIs(True)

    checkPPIGen = map(wrapperCheckPPIs, subsampleSizes,  
        repeats, repeat(proteinsData), repeat(corum))  # First For Cycle

    allAUC = dict(checkPPIGen)

    # DEPRECATED UNPARARELIZED CODE
    # # Random Sampling
    # for sampleNum in subsampleSizes:

    #     iteration = 0
    #     aucList = list()
    #     sumOfTime = 0
    #     while iteration < repeats:

    #         proteinsDatSample = proteinsData.sample(n = sampleNum, axis=0, random_state=RANDOMSTATE)
    #         pairwiseCorr = utils.getPairwiseCorrelation(
    #             proteinsDatSample, None, str(sampleNum) + ' samples', False)
    #         del proteinsDatSample

    #         pairwiseCorr = utils.addGroundTruth(corum, pairwiseCorr, 'corum', None)
    #         end = time.time()
    #         sumOfTime += end-start

    #         corrCumSum = np.cumsum(pairwiseCorr['corum']) / np.sum(pairwiseCorr['corum'])
    #         indexes = np.array(pairwiseCorr.reset_index().index) / pairwiseCorr.shape[0]
    #         AUC = auc(indexes, corrCumSum)
    #         aucList.append(AUC)
    #         del pairwiseCorr
    #         iteration += 1

    #     print(aucList)
    #     print(sumOfTime/repeats)
    #     allAUC[str(sampleNum)] = aucList

    # Plot each AUC in the various boxplot chart
    corumAUC = 0.76
    fig, ax = plt.subplots(figsize=(40, 8))
    ax.boxplot(allAUC.values(), labels=allAUC.keys())
    ax.set_ybound(lower=0.2, upper=1)
    ax.set_ylabel("AUC", fontsize=14)
    ax.set_xlabel("Sampling Number", fontsize=14)
    ax.axhline(y=corumAUC, color='red', linestyle='-', label='BaseModel AUC')
    ax.axhline(y=0.9*corumAUC, color='blue', linestyle=':', label='90% of BaseModel AUC')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)


    plt.savefig(PATH  +  "/images/aucPerSamplingNumberv2.0.png",
                bbox_inches="tight")


if __name__ == '__main__':
    subsamplingList = list(range(5,950,5))
    repeatsList= [round(900/repeat) + 5 if round(900/repeat) >= 4 and round(900/repeat) <= 100 else 100 if round(900/repeat)*2 > 100 else 5  for repeat in subsamplingList]
    start = time.time()
    randomSubSamplingAUC(proteinsData, subsamplingList, repeatsList)
    end = time.time()
    sumOfTime = end-start
    print(sumOfTime)
