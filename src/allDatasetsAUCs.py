import pandas as pd
import time as t
from resources import *





if __name__ == "__main__":


    ogAndVAEModel: PairwiseCorrMatrix = read(PATH + '/datasetsTese/ogAndVAEPairwiseCorrs.pickle.gz')
    ogAndVAEData: pd.DataFrame = ogAndVAEModel.data


    corum:ppiDataset = read(PATH + '/externalDatasets/corum.pickle.gz')
    biogrid:ppiDataset = read(PATH + '/externalDatasets/biogrid.pickle.gz')
    string150:ppiDataset = read(PATH + '/externalDatasets/string150.pickle.gz')
    string400:ppiDataset = read(PATH + '/externalDatasets/string400.pickle.gz')
    string700:ppiDataset = read(PATH + '/externalDatasets/string700.pickle.gz')
    string900:ppiDataset = read(PATH + '/externalDatasets/string900.pickle.gz')

    ppisSets:list[set[str,str]] = [corum.ppis, biogrid.ppis, string150.ppis, string400.ppis, string700.ppis, string900.ppis]
    datasetNames:list[str] = ['corum', 'biogrid', 'string Low', 'string Moderate', 'string High', 'string Highest']


    for ppis, datasetName in zip(ppisSets, datasetNames):

        start = t.time()   
        ogAndVAEModel.addGroundTruthNeo(ppis, datasetName)
        print(t.time() - start)
    ogAndVAEModel.write(PATH + '/datasetsTese/ogAndVAEPairwiseCorrs.pickle.gz')
