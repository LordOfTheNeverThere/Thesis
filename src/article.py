    # Imports
    from typing import Iterable
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from resources import read, PairwiseCorrMatrix ,ProteinsMatrix, PATH, drawRecallCurves
    import numpy as np
    from sklearn.linear_model import LinearRegression
    from scipy.stats import shapiro, pearsonr, chi2_contingency
    from statsmodels.stats.diagnostic import het_white, het_breuschpagan


if __name__ == '__main__':

    #To create Images for Article

    proteomics  :ProteinsMatrix= read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    vaeProteomics :ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')
    meanProteomicsLib :ProteinsMatrix = read(PATH + '/internal/proteomics/mean75PVProteomics.pickle.gz')
    meanProteomicsConser :ProteinsMatrix = read(PATH + '/internal/proteomics/mean80PVProteomics.pickle.gz')
    instances = [vaeProteomics, meanProteomicsLib, meanProteomicsConser]

    for inst in instances:

        inst.shapiroWilksTest()
        whitened ,_ =ProteinsMatrix.whitening(inst, inst, True)
        whitened = ProteinsMatrix(None, whitened)
        whitened.shapiroWilksTest()



    for inst in instances:
        inst.independenceSamples()