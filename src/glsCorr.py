import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset, PairwiseCorrMatrix


if __name__ == '__main__':
        
    vaeProteomics = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')

    numPC:int = 10
    factorsName:str=''

