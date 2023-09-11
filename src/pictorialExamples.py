import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset, PairwiseCorrMatrix


if __name__ == '__main__':


    def getExample(
            proteinsCorrelated:bool, 
            normallyHighDR:bool, 
            clusterAbove:bool, 
            lineSamplesNum:int = 100, 
            clusterSamplesNum:int = 20):
        np.random.seed(0)
        
        results ={'index':[], 'Px':[], 'Py':[], 'DR':[]}
        lineSamples = range(lineSamplesNum)
        clusterSamples = range(clusterSamplesNum)

        if proteinsCorrelated:
            Px = np.random.normal(0, 1, lineSamplesNum)
            Px = Px + np.random.normal(0, 0.1, lineSamplesNum)
            Py = 1.5*Px + np.random.normal(0, 0.1, lineSamplesNum)

            PxCluster = np.random.normal(4, 1, clusterSamplesNum)
            PxCluster = PxCluster + np.random.normal(0, 0.5, clusterSamplesNum)

            if clusterAbove:
                
                PyCluster = 1.5*PxCluster + np.random.normal(0, 0.5, clusterSamplesNum)

            else:

                PyCluster = 1.5*PxCluster + np.random.normal(0, 0.5, clusterSamplesNum) - 7
        else:
            Px = np.random.normal(0, 1, lineSamplesNum)
            Px = Px + np.random.normal(0, 0.1, lineSamplesNum)
            Py = -1.5*Px - np.random.normal(0, 0.1, lineSamplesNum)

            PxCluster = np.random.normal(4, 1, clusterSamplesNum)
            PxCluster = PxCluster + np.random.normal(0, 0.5, clusterSamplesNum)

            if clusterAbove:
                
                PyCluster = -1.5*PxCluster - np.random.normal(0, 0.5, clusterSamplesNum)

            else:

                PyCluster = -1.5*PxCluster - np.random.normal(0, 0.5, clusterSamplesNum) - 7

        Px = np.concatenate((Px, PxCluster))
        Py = np.concatenate((Py, PyCluster))
        
        if normallyHighDR:
            DR = np.random.normal(10, 1, lineSamplesNum); DRCluster =  np.random.normal(5, 1, clusterSamplesNum)
        else:
            DR = np.random.normal(5, 1, lineSamplesNum); DRCluster =  np.random.normal(10, 1, clusterSamplesNum)

        DR = np.concatenate((DR, DRCluster))

        results['index'] = list(range(lineSamplesNum + clusterSamplesNum))
        results['Px'] = [Px]
        results['Py'] = [Py]
        results['DR'] = [DR]


        return results

        
# Sample 100 samples from the normal distribution, with mean 0 and std 1, adding some noise
    #1. One where the low IC50 is the linear line and the high IC50 are cluster above

    getExample(proteinsCorrelated=True, normallyHighDR=False, clusterAbove=True, clusterSamplesNum=5, lineSamplesNum=15)
    

