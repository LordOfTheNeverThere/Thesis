# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix, ppiDataset
import pickle
import gzip
import utils

from env import PATH






if __name__ == '__main__':
    ogPairwise = utils.read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')
    glsPairwise = utils.read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    utils.drawRecallCurves([ogPairwise,glsPairwise], ['blue', 'red'],'ogVsGLSRecallCurve.png')
    # 