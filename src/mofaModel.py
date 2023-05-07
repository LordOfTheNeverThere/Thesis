# %% Imports
import pandas as pd
import seaborn as sb
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import pickle
import gzip
from resources import *

#Load Dataset
def mofaBaseModel():
    with gzip.open(PATH + '/datasetsTese/mofaPairwiseCorr.pickle.gz', 'rb') as file:
        mofaPairwise = pickle.load(file)
    file.close()

    with gzip.open(PATH + '/datasetsTese/baseModelFiltered.pickle.gz', 'rb') as file:
        ogPairwise = pickle.load(file)
    file.close()


    utils.drawRecallCurves([mofaPairwise, ogPairwise], [
                           'blue', 'red'], "mofaVsOgRecallCurve.png")


def mofa2():
    """In this model we are only using the proteins in the mofa proteomics matrix which are also on the og matrix"""
    with gzip.open(PATH + '/datasetsTese/mofaPairwiseCorr.pickle.gz', 'rb') as file:
        mofaPairwise = pickle.load(file)
    file.close()

    with gzip.open(PATH + '/datasetsTese/baseModelFiltered.pickle.gz', 'rb') as file:
        ogPairwise = pickle.load(file)
    file.close()
    ogPairwise.data = ogPairwise.data.dropna()
    indexesOfInterest = mofaPairwise.data.index.intersection(ogPairwise.data.index) 
    mofaPairwise.data = mofaPairwise.data.loc[indexesOfInterest]

    mofaPairwise.aucCalculator('corum', 'MofaWithOverlap')


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
            pairwiseCorr.addGroundTruth(corum, 'corum', None)
            pairwiseCorr.aucCalculator('corum')

            allPairwiseCorrs.append(pairwiseCorr)



        
        hues = ['red', 'blue', 'green', 'black', 'purple', 'brown']
        utils.drawRecallCurves(allPairwiseCorrs, hues, PATH + "/images/mofaRecallPresenceThreshv3.png")




def opposingIntensities():

    with gzip.open(PATH + '/datasetsTese/mofaPairwiseCorr.pickle.gz', 'rb') as f:
        mofa = pickle.load(f)
    f.close()

    with gzip.open(PATH + '/datasetsTese/baseModelFiltered.pickle.gz', 'rb') as f:
        baseModel = pickle.load(f)
    f.close()

    with gzip.open(PATH + '/datasetsTese/mofaProteomics.pickle.gz', 'rb') as f:
        mofaProteomics = pickle.load(f)
    f.close()

    

    queriedFrame = baseModel.compare(mofa, 'corum == 1','corum == 1')
    queriedFrame['corrDiference'] = abs(queriedFrame['globalCorrelation'] - queriedFrame['mofaCorrelation'])
    queriedFrame.sort_values(by='corrDiference', ascending=False, inplace=True)
    highestDiference = queriedFrame.head(5)
    corrDiference = list(highestDiference['corrDiference'])
    indexes = list(highestDiference.index)
    setOfPPIs = [(proteins.split(';')[0], proteins.split(';')[1] ) for proteins in indexes]#Unpack PPIs of opposing inensities into tuples of proteins 

    fig, ax = plt.subplots(5, 1, figsize=(15, 40))
    for index, ppi in enumerate(setOfPPIs):
        proteinA = pd.Series(mofaProteomics.data[ppi[0]])
        proteinB = pd.Series(mofaProteomics.data[ppi[1]])

        ax[index].scatter(proteinA, proteinB)
        ax[index].set_xlabel(ppi[0])
        ax[index].set_ylabel(ppi[1])
        ax[index].set_title('Δcorr ==' + str(corrDiference[index]) + '-' + ppi[0] + ';' +ppi[1])
        ax[index].tick_params(labelsize=16)

    plt.savefig('../images/pxVSpyHighestDifference.png')




 

if __name__ == '__main__':


    opposingIntensities()
