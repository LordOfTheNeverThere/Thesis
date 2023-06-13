import pandas as pd
import time as t
from resources import *





if __name__ == "__main__":


    ogAndVAEModel: PairwiseCorrMatrix = read(PATH + '/internal/ogAndVAEPairwiseCorrs.pickle.gz')
    ogAndVAEData: pd.DataFrame = ogAndVAEModel.data

    ogAndVAEData.to_csv(PATH + '/internal/ogAndVAEPairwiseCorrs.csv.gz', index = False, compression='gzip')


    corum:ppiDataset = read(PATH + '/external/corum.pickle.gz')
    biogrid:ppiDataset = read(PATH + '/external/biogrid.pickle.gz')
    string150:ppiDataset = read(PATH + '/external/string150.pickle.gz')
    string400:ppiDataset = read(PATH + '/external/string400.pickle.gz')
    string700:ppiDataset = read(PATH + '/external/string700.pickle.gz')
    string900:ppiDataset = read(PATH + '/external/string900.pickle.gz')

    ppisSets:list[set[str,str]] = [corum.ppis, biogrid.ppis, string150.ppis, string400.ppis, string700.ppis, string900.ppis]
    datasetNames:list[str] = ['corum', 'biogrid', 'string Low', 'string Moderate', 'string High', 'string Highest']


    for ppis, datasetName in zip(ppisSets, datasetNames):

        start = t.time()   
        ogAndVAEModel.addGroundTruth(ppis, datasetName)
        print(t.time() - start)
    ogAndVAEModel.write(PATH + '/internal/ogAndVAEPairwiseCorrs.pickle.gz')
