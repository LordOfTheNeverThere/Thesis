import pandas as pd
from resources import PATH, CPUS, read, DrugResponseMatrix, GeneDependency, ResiduesMatrix, ProteinsMatrix, PairwiseCorrMatrix
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

    #Proteomics

    ogProteomics = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
    vaeProteomics = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz')
    pv75MeanProteomics = read(PATH + '/internal/proteomics/mean75%PVProteomics.pickle.gz')
    pv80MeanProteomics = read(PATH + '/internal/proteomics/mean80%PVProteomics.pickle.gz')
    #How was pv75MeanProteomics created?
    pv75MeanDf:pd.DataFrame = ogProteomics.data.copy()
    pv75MeanDf = pv75MeanDf.dropna(axis=1, thresh=round(0.75*pv75MeanDf.shape[0]))
    pv75MeanDf = pv75MeanDf.fillna(dict(pv75MeanDf.mean()))



    #Drug Response
    drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')


    #Samplesheet
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)

    #Gene Dependency
    geneDependency:GeneDependency = read(PATH + '/internal/geneInteractionModel/geneDependency.pickle.gz')

    


