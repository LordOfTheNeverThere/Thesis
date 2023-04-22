# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix, ppiDataset
import matplotlib.pyplot as plt

from env import PATH




if __name__ == '__main__':

    # Pairwise Corr
    pairwiseCorrPValues = PairwiseCorrMatrix(PATH + '/datasetsTese/baseModelPairwiseWithpValues.csv.gz', compression = 'gzip', index_col='PPI')
    # External ppi data
    corum = ppiDataset(PATH + '/externalDatasets/corumPPI.csv.gz', compression='gzip')
    biogrid = ppiDataset(PATH + '/externalDatasets/biogridPPIHuman.csv.gz', compression='gzip')
    string = ppiDataset(PATH + '/externalDatasets/stringPPI900Selected.csv.gz', compression='gzip')
    
    #Add ground truths to the pairwise coor matrix
    for dataset in [(corum,'corum'), (biogrid,'biogrid'), (string,'string')]:
        pairwiseCorrPValues.addGroundTruth(dataset[0].getPPIs(dataset[1]), dataset[1], PATH + '/datasetsTese/baseModePairwiseWithpValues.csv.gz' if dataset[0] == 'string' else None)
    
