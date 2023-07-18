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
    vaeGLSPairwise.data['fdr'] = multipletests(vaeGLSPairwise.data['p-value'], method='fdr_bh')[1]
    ppisOfInterest = set(vaeGLSPairwise.data.query("corum ==1 and fdr < 0.01").index)
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
    dummy.filepath = PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/PCARegressor.pickle.gz'
    dummy.write()


    # dummy:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/fdrPerPPIRegressor.pickle.gz')        
    # data = dummy.data
    # data = data.loc[data['info']['fdr'] < 0.01]
    # betaThresh = data['effectSize']['interaction'].quantile(0.98) # define a beta threshold based on a quantile given by the use
    # data = data.loc[abs(data['effectSize']['interaction']) > betaThresh] # subset data to only include betas above the threshold
    # data = data.sort_values(by=[('info','logLikePValue')], ascending=[True])

    # dummy.volcanoPlot('volcanoPlotDrInteractionPxModelFDRBuffer.png') # 3579956 points
    # drugRes.data = drugRes.data.T
    # dummy.scatterTheTopVolcano('topVolcanoPlotScatter.png', ogProteomics, drugRes, topNumber=10)

