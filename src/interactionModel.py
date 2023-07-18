import pandas as pd
from resources import PATH, CPUS, read, DrugResponseMatrix, DRInteractionPxModel, ResiduesMatrix, ProteinsMatrix, PairwiseCorrMatrix
import statsmodels.api as sm
import statsmodels.formula.api as smf
from typing import Iterable
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.multitest import multipletests
import time as t

        

if __name__ == '__main__':

    #py ~ M + px + drugres + px:drugres

    drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    drugRes.data = drugRes.data.T
    samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    vaeProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz') #used for PCA computation
    ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz') #used for the interaction model class
           
    # using only corum ppis that we were able to recall, with high confidence
    vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')
    ppisOfInterest = set(vaeGLSPairwise.data.query("corum ==1 and fdr < 0.01").head(3).index)
    ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}

    #Cofounding Factors, use The samplesheet's growth properties or the 10 PC of the vaeProteomics dataframe    
    growthProps = pd.get_dummies(samplesheet['growth_properties'])
    growthProps = growthProps.rename(columns={'Semi-Adherent': 'SemiAdherent'})

    pcFactors = vaeProteomics.PCA(filepath='pca.png', factorsName='PC').factors

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)


    dummy = DRInteractionPxModel(ppisOfInterest, ogProteomics, drugRes, pcFactors)
    fit = dummy.fit()
    # dummy.filepath = PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/fdrPerPPIRegressor.pickle.gz'
    # dummy.write()


    dummy:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/fdrPerPPIRegressor.pickle.gz')    
    dummy.volcanoPlot('volcanoPlotDrInteractionPxModel.png') # 142393 points
    drugRes.data = drugRes.data.T
    dummy.scatterTheTopVolcano('topVolcanoPlotScatter.png', ogProteomics, drugRes)

