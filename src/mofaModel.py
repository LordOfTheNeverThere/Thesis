# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import utils
from classes import ppiDataset, ProteinsMatrix, PairwiseCorrMatrix

from env import PATH

#Load Dataset
def mofaBaseModel():
    # proteinsData = ProteinsMatrix(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', index_col='Unnamed: 0')

    # #  Load external Datasets

    # corum = ppiDataset(filepath=PATH + '/externalDatasets/corumPPI.csv.gz')
    # corum = corum.getPPIs('corum')

    # pairwiseCorr = proteinsData.pearsonCorrelations(None, columnName='mofaCorrelation', counting=False)
    pairwiseCorr = PairwiseCorrMatrix(PATH + '/datasetsTese/baseMOFAPairwiseCorr.csv.gz', data=None, compression='gzip', index_col='PPI')
    ogPairwise = PairwiseCorrMatrix( PATH+'/datasetsTese/baseModel.csv.gz', data=None, compression='gzip', index_col='PPI')
    print(pairwiseCorr.data)
    print(ogPairwise.data)

    # pairwiseCorr.addGroundTruth(corum, 'corum', PATH + '/datasetsTese/baseMOFAPairwiseCorr')

    pairwiseCorr.aucCalculator('corum', 'mofaModel')
    ogPairwise.aucCalculator('Corum', 'baseModel')

    utils.drawRecallCurves([pairwiseCorr, ogPairwise], ['blue', 'red'],"../images/mofaVsOgRecallCurve.png")


def mofa2():
    """In this model we are only using the proteins in the mofa proteomics matrix which are also on the og matrix"""
    mofaPairwise = PairwiseCorrMatrix(PATH + '/datasetsTese/baseMOFAPairwiseCorr.csv.gz', data=None, compression='gzip', index_col='PPI')
    ogPairwise = PairwiseCorrMatrix(PATH+'/datasetsTese/baseModel.csv.gz', data=None ,compression='gzip',index_col='PPI')
    ogPairwise.data = ogPairwise.data.dropna()
    indexesOfInterest = mofaPairwise.data.index.intersection(ogPairwise.data.index) 
    mofaPairwise.data = mofaPairwise.data.loc[indexesOfInterest]
    print(mofaPairwise.data)

    mofaPairwise.aucCalculator('Corum', 'MofaWithOverlap')


    utils.drawRecallCurves([mofaPairwise], ['red'], 'reCurveMofaOgFiltered.png')



    
def mofa3(threshold: float | list[float]):

    """Similarly to the above function, we will filter proteins of the mofa matrix with those of the og matrix which have at least a threshold of presistence
    Args:
        threshold (float): Filtering threshold which states the value of presence a protein must have in order for it to be included in the final mofa matrix.
        If we want proteins that are present in at leat 10% of rows we set it to 0.10
    """
    mofaProteins = ProteinsMatrix(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', {'index_col':'Unnamed: 0'})
    ogProteins = ProteinsMatrix(PATH+'/datasetsTese/proteomicsDataTrans.csv', {'index_col': 'modelID'})
    
    assert type(threshold) == float or type( threshold) == list, 'Wrong Type for prop threshold'

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
        utils.drawRecallCurves(allPairwiseCorrs, hues, PATH + "/images/mofaRecallPresenceThreshv3.png")




def opposingIntensities():
    mofa = PairwiseCorrMatrix(
        PATH + '/datasetsTese/baseMOFAPairwiseCorr.csv.gz', compression='gzip')
    baseModel = PairwiseCorrMatrix(PATH + '/datasetsTese/BaseModel.csv.gz', compression='gzip')
    print(baseModel.data)
    print(mofa.data)

    queriedFrame = baseModel.compare(mofa, 'Corum == 1','corum == 1')

    print(queriedFrame)
    queriedFrame['corrDiference'] = abs(queriedFrame['globalCorrelation'] - queriedFrame['mofaCorrelation'])
    highestDiference = queriedFrame.sort_values(by='corrDiference', ascending=False, inplace=True).head(5)
    indexes = list(highestDiference.index)
    setOfPPIs = {(proteins.split(';')[0], proteins.split(';')[1]) for proteins in indexes} #Unpack PPIs of opposing inensities into tuples of proteins

    for ppi in setOfPPIs:
        pass



 

if __name__ == '__main__':
    mofaBaseModel()
    #opposingIntensities()

