# Imports
import pandas as pd
import matplotlib.pyplot as plt
from resources import *







if __name__ == '__main__':

    vaePearPair:PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')
    vaePearPair.aucsCalculator(['corum', 'corum'], 'VAE Pearson', ['pearsonR', 'pValue'], [False, True], PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')