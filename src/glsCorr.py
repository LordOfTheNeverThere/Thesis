import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt
import time as t
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset


if __name__ == '__main__':
    vaeProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')
