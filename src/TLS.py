
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from resources import *



proteomics:ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
glsPairwise:PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
print(glsPairwise)

proteomics.tlsResidues(glsPairwise) # Correct This!


