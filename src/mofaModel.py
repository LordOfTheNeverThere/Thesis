# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from classes import ppiDataset, ProteinsMatrix, PairwiseCorrMatrix

from env import PATH

# %% Load Dataset
def mofaBaseModel():
    proteinsData = ProteinsMatrix(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', index_col='Unnamed: 0')

    #  Load external Datasets

    corum = ppiDataset(filepath=PATH + '/externalDatasets/corumPPI.csv.gz')
    corum = corum.getPPIs(True)

    pairwiseCorr = proteinsData.pearsonCorrelations(
        filepath=None, columnName='mofaCorrelation', counting=False)
    
    pairwiseCorr = PairwiseCorrMatrix(data = pairwiseCorr)

    pairwiseCorrData = pairwiseCorr.addGroundTruth(corum, 'Corum', None)

    pairwiseCorrData.to_csv(PATH + '/datasetsTese/baseMOFACorr.csv.gz', compression='gzip')

    # pairwiseCorr.aucCalculator('Corum')

    # utils.drawRecallCurves([pairwiseCorr], ['blue'],PATH + "/images/mofaModelRecallCurve.png")


# def mofa2(): DDEPRECATED Filtering did nothing
#     """In this model we are only using the proteins in the mofa proteomics matrix which are also on the og matrix
#     """
#     mofaProteins = pd.read_csv(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', compression='gzip', index_col='Unnamed: 0')
#     ogProteins = pd.read_csv(PATH+'/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')
    
#     # Filter DataFrame based on columns from the list, even if the list contains values not present in the DataFrame
#     filteredMofaProteins = mofaProteins[mofaProteins.columns.intersection(ogProteins.columns)]
#     #The matrix suffers now change soo this is rather usless, keep it if however you wish to do some future tweaking
    
def mofa3(threshold: float | list[float]):

    """Similarly to the above function, we will filter proteins of the mofa matrix with those of the og matrix which have at least a threshold of presistence
    Args:
        threshold (float): Filtering threshold which states the value of presence a protein must have in order for it to be included in the final mofa matrix.
        If we want proteins that are present in at leat 10% of rows we set it to 0.10
    """
    mofaProteins = ProteinsMatrix(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', {'index_col':'Unnamed: 0'})
    ogProteins = ProteinsMatrix(PATH+'/datasetsTese/proteomicsDataTrans.csv', {'index_col': 'modelID'})
    

    if type(threshold) == float:
        threshold = list(threshold)
    if type(threshold) == list:

        allPairwiseCorrs = list()

        for thresh in threshold:
            ogProteins.data.dropna(axis=1, thresh=round(ogProteins.shape[0] * thresh), inplace=True) 
            # We are dropping columns with at least a threshold of missingness, hence we are subjugating the matrix to mandatory threshold of presence, confusing... I know
            mofaProteins.data = mofaProteins.data[mofaProteins.columns.intersection(ogProteins.columns)]
            #  Load external Datasets

            corum = ppiDataset(filepath=PATH + '/externalDatasets/corumPPI.csv.gz')
            corum = corum.getPPIs(True)

            pairwiseCorr = mofaProteins.pearsonCorrelations(filepath=None, columnName='mofaCorrelation', counting=False)
            pairwiseCorr = PairwiseCorrMatrix(data=pairwiseCorr)
            pairwiseCorr.addGroundTruth(corum, 'Corum', None)
            pairwiseCorr.aucCalculator('Corum')

            allPairwiseCorrs.append(pairwiseCorr)



        
        hues = ['red', 'blue', 'green', 'black', 'purple', 'brown']
        utils.drawRecallCurves(allPairwiseCorrs, hues,
                               PATH + "/images/mofaRecallPresenceThreshv3.png")

        
    else:
        print('Wrong Type for prop threshold')



def opposingIntensities():
    mofa = PairwiseCorrMatrix(PATH + '/datasetsTese/baseMOFACorr.csv.gz')
    baseModel = PairwiseCorrMatrix(PATH + '/datasetsTese/BaseModelPairwise.csv', gziped=False)

    queriedFrame = baseModel.compare(mofa, 'Corum == 1 & globalCorrelation  > 0.8','Corum == 1 & mofaCorrelation  < 0.4')

    print(queriedFrame)
 

if __name__ == '__main__':
    # mofa3([0.45, 0.46, 0.47, 0.48, 0.49, 0.5])
    mofaBaseModel()
    # opposingIntensities()
