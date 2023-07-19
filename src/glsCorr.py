import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt

from statsmodels.stats.multitest import multipletests
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset, PairwiseCorrMatrix


if __name__ == '__main__':

    glsVae : PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')
    glsVae.data.rename(columns={'pValue': 'p-value', 'beta':'coef'}, inplace=True)
    glsVae.data['fdr'] = multipletests(glsVae.data['p-value'], method='fdr_bh')[1]

    pearsonVae : PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/pearsonPairCorr.pickle.gz')
    pearsonVae.data.rename(columns={'pValue': 'p-value', 'pearsonR':'coef'}, inplace=True)
    pearsonVae.data['fdr'] = multipletests(pearsonVae.data['p-value'], method='fdr_bh')[1]
    Done#

    glsMean75: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')
    glsMean75.data['fdr'] = multipletests(glsMean75.data['p-value'], method='fdr_bh')[1]

    pearsonMean75: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr75PV.pickle.gz')
    pearsonMean75.data.rename(columns={'pValue': 'p-value', 'pearsonR':'coef'}, inplace=True)
    pearsonMean75.data['fdr'] = multipletests(pearsonMean75.data['p-value'], method='fdr_bh')[1]


    glsMean80: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr80PV.pickle.gz')
    glsMean80.data.rename(columns={'pValue': 'p-value', 'beta':'coef'}, inplace=True)
    glsMean80.data['fdr'] = multipletests(glsMean80.data['p-value'], method='fdr_bh')[1]

    pearsonMean80: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr80PV.pickle.gz')
    pearsonMean80.data.rename(columns={'pValue': 'p-value', 'pearsonR':'coef'}, inplace=True)
    pearsonMean80.data['fdr'] = multipletests(pearsonMean80.data['p-value'], method='fdr_bh')[1]

    pearsonOg: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/OG/baseModelFiltered.pickle.gz')
    pearsonOg.data.rename(columns={'pValue': 'p-value', 'pearsonR':'coef'}, inplace=True)
    pearsonOg.data['fdr'] = multipletests(pearsonOg.data['p-value'], method='fdr_bh')[1]
    
    pairwiseList = [glsMean75, pearsonMean75, glsMean80, pearsonMean80, pearsonOg]
    
    #PairwiseCorrMatrix.addGroundTruths(pairwiseList)

    PairwiseCorrMatrix.getAucs(pairwiseList)