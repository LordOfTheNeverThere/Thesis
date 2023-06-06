# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *
import time as t







if __name__ == '__main__':

    vaePearPair:PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    start = t.time()
    ppis = read(PATH + '/externalDatasets/corum.pickle.gz').ppis
    vaePearPair.addGroundTruth(ppis, 'corum')
    print(t.time() - start)
    start = t.time()
    vaePearPair.aucsCalculator(['corum', 'corum'], 'VAE Pearson', ['pearsonR', 'pValue'], [False, True], PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    print(t.time() - start)  # print elapsed time