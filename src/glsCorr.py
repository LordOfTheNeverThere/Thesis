import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt

from statsmodels.stats.multitest import multipletests
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset, PairwiseCorrMatrix


if __name__ == '__main__':

    glsVae : PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')

    pearsonVae : PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/pearsonPairCorr.pickle.gz')

    glsMean75: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')

    pearsonMean75: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr75PV.pickle.gz')

    glsMean80: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr80PV.pickle.gz')

    pearsonMean80: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr80PV.pickle.gz')

    pearsonOg: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/OG/baseModelFiltered.pickle.gz')
    
    pairwiseList = [glsVae, pearsonVae, glsMean75, pearsonMean75, glsMean80, pearsonMean80, pearsonOg]
    
    #PairwiseCorrMatrix.addGroundTruths(pairwiseList)

    PairwiseCorrMatrix.getAucs(pairwiseList)