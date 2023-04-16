# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from classes import ppiDataset

from env import PATH

# %% Load Dataset
def mofaBaseModel():
    proteinsData = pd.read_csv(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', compression='gzip', index_col='Unnamed: 0')

    #  Load external Datasets

    corum = ppiDataset(filename=PATH + '/externalDatasets/corumPPI.csv.gz')
    corum = corum.getPPIs(True)

    pairwiseCorr = utils.getPairwiseCorrelation(
        fileName=None, data=proteinsData, columnName='mofaCorrelation')

    pairwiseCorr = utils.addGroundTruth(corum, pairwiseCorr, 'Corum', None)

    pairwiseCorr.to_csv(PATH + '/datasetsTese/baseMOFACorr.csv.gz', compression='gzip', index= False)

    # Create Recall Curves
    # corrCumSum = np.cumsum(
    #     pairwiseCorr['Corum']) / np.sum(pairwiseCorr['Corum'])
    # indexes = np.array(pairwiseCorr.reset_index().index) / \
    #     pairwiseCorr.shape[0]
    # AUC = auc(indexes, corrCumSum)


    # _, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)
    # ax.plot(
    #     indexes,
    #     corrCumSum,
    #     label=f"(AUC {AUC:.2f})",
    #     c='blue',
    # )

    # ax.plot([0, 1], [0, 1], "k--", lw=0.3)
    # ax.legend(loc="lower right", frameon=False)

    # ax.set_ylabel("Cumulative sum")
    # ax.set_xlabel("Ranked correlation")
    # ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

    # # plt.savefig(f"{RPATH}/PPInteractions_roc_curves_overlap.pdf",
    # #             bbox_inches="tight")
    # plt.savefig("mofaModelRecallCurve.png",
    #             bbox_inches="tight")
    # plt.close("all")

# def mofa2(): DDEPRECATED Filtering did nothing
#     """In this model we are only using the proteins in the mofa proteomics matrix which are also on the og matrix
#     """
#     mofaProteins = pd.read_csv(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', compression='gzip', index_col='Unnamed: 0')
#     ogProteins = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')
    
#     # Filter DataFrame based on columns from the list, even if the list contains values not present in the DataFrame
#     filteredMofaProteins = mofaProteins[mofaProteins.columns.intersection(ogProteins.columns)]
#     #The matrix suffers now change soo this is rather usless, keep it if however you wish to do some future tweaking
    
def mofa3(threshold: float or list[float]):

    """Similarly to the above function, we will filter proteins of the mofa matrix with those of the og matrix which have at least a threshold of presistence
    Args:
        threshold (float): Filtering threshold which states the value of presence a protein must have in order for it to be included in the final mofa matrix.
        If we want proteins that are present in at leat 10% of rows we set it to 0.10
    """

    mofaProteins = pd.read_csv(
        PATH + '/datasetsTese/proteomicsMOFA.csv.gz', compression='gzip', index_col='Unnamed: 0')
    ogProteins = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')
    

    if type(threshold) == float:
        threshold = list(threshold)
    if type(threshold) == list:

        allAUC = dict()

        for thresh in threshold:
            ogProteins.dropna(axis=1, thresh=round(ogProteins.shape[0] * thresh), inplace=True) 
            # We are dropping columns with at least a threshold of missingness, hence we are subjugating the matrix to mandatory threshold of presence, confusing... I know
            filteredMofaProteins = mofaProteins[mofaProteins.columns.intersection(ogProteins.columns)]
            #  Load external Datasets

            corum = ppiDataset(filename=PATH + '/externalDatasets/corumPPI.csv.gz')
            corum = corum.getPPIs(True)

            pairwiseCorr = utils.getPairwiseCorrelation(
                fileName=None, data=filteredMofaProteins, columnName='mofaCorrelation')

            pairwiseCorr = utils.addGroundTruth(corum, pairwiseCorr, 'Corum', None)

            # Create Recall Curves
            corrCumSum = np.cumsum(
                pairwiseCorr['Corum']) / np.sum(pairwiseCorr['Corum'])
            indexes = np.array(pairwiseCorr.reset_index().index) / \
                pairwiseCorr.shape[0]
            AUC = auc(indexes, corrCumSum)
            allAUC[str(thresh)] = {'corrCumSum': corrCumSum,'indexes': indexes,'auc': AUC}


        _, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)
        hues = ['red', 'blue', 'green', 'black', 'purple', 'brown']

        for index,key in enumerate(allAUC):

            ax.plot(
                allAUC[key]['indexes'],
                allAUC[key]['corrCumSum'],
                label=f"(AUC {allAUC[key]['auc']:.2f}) thresh=" + key,
                c=hues[index],
            )

        ax.plot([0, 1], [0, 1], "k--", lw=0.3)
        ax.legend(loc="lower right", frameon=False)

        ax.set_ylabel("Cumulative sum")
        ax.set_xlabel("Ranked correlation")
        ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

        # plt.savefig(f"{RPATH}/PPInteractions_roc_curves_overlap.pdf",
        #             bbox_inches="tight")
        plt.savefig("mofaRecallPresenceThreshv3.png",
                    bbox_inches="tight")
        plt.close("all")

        
    else:
        print('Wrong Type for prop threshold')

if __name__ == '__main__':
    # mofa3([0.45, 0.46, 0.47, 0.48, 0.49, 0.5])
    mofaBaseModel()
