# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix, ppiDataset
import matplotlib.pyplot as plt
import multiprocessing as mp

from env import PATH

CPUS = 6
assert CPUS < mp.cpu_count() - 1




if __name__ == '__main__':

    # Pairwise Corr
    pairwiseCorrPValues = PairwiseCorrMatrix(PATH + '/datasetsTese/baseModel.csv.gz', compression = 'gzip', index_col='PPI')
    # External ppi data
    corum = ppiDataset(PATH + '/externalDatasets/corumPPI.csv.gz', compression='gzip')
    biogrid = ppiDataset(PATH + '/externalDatasets/biogridPPIHuman.csv.gz', compression='gzip')
    string = ppiDataset(PATH + '/externalDatasets/stringPPI900Selected.csv.gz', compression='gzip')

    def addGroundTruthWrapper(dataset: tuple):
        print("dataset == " + dataset[1])
        pairwiseCorrPValues.addGroundTruth(dataset[0].getPPIs(
            dataset[1]), dataset[1], None)
    
    datasets = [(corum, 'corum'), (biogrid, 'biogrid'), (string, 'string')]
    # Add ground truths to the pairwise coor matrix
    with mp.Pool(CPUS) as process:
        process.map(addGroundTruthWrapper, datasets)

    pairwiseCorrPValues.data.to_csv(
        PATH + '/datasetsTese/baseModel.csv.gz', compression='gzip')
        
    
