import pandas as pd
from resources import PATH, CPUS, read, DrugResponseMatrix, DRInteractionPxModel, ResiduesMatrix, ProteinsMatrix, PairwiseCorrMatrix
import statsmodels.api as sm
import statsmodels.formula.api as smf
from typing import Iterable
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.multitest import multipletests
import time as t

# main
if __name__ == '__main__':

    # Pairwise Corrs
    #Mean
    pv75Mean = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')
    pv80Mean = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr80PV.pickle.gz')