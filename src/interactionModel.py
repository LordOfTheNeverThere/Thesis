# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as t
from statsmodels.stats.multitest import multipletests
from resources import GeneDependency, DRPxPyInteractionPxModel, ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix


        

if __name__ == '__main__':

    #py ~ M + px + drugres + px:drugres

    # drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
    # drugRes.data = drugRes.data.T
    # samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
    # vaeProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAE.pickle.gz') #used for PCA computation
    # ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz') #used for the interaction model class
           
    # # using only corum ppis that we were able to recall, with high confidence
    # vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')
    # vaeGLSPairwise.data['fdr'] = multipletests(vaeGLSPairwise.data['p-value'], method='fdr_bh')[1]
    # ppisOfInterest = set(vaeGLSPairwise.data.query("corum ==1 and fdr < 0.01").index)
    # ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}

    # #Cofounding Factors, use The samplesheet's growth properties or the 10 PC of the vaeProteomics dataframe    
    # growthProps = pd.get_dummies(samplesheet['growth_properties'])
    # growthProps = growthProps.rename(columns={'Semi-Adherent': 'SemiAdherent'})

    # pca, pcFactors = vaeProteomics.PCA(factorsName='PC', numPC=5)

    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.width', None)
    # pd.set_option('display.max_colwidth', -1)

    # dummy = DRPxPyInteractionPxModel(ppisOfInterest, ogProteomics, drugRes.data, pcFactors)
    # start = t.time()
    # fit = dummy.fit(numOfCores = 38)
    # dummy.filepath = PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/neoDrugRegressor.pickle.gz'
    # dummy.write()
    # print(f'fitting took {t.time() - start} seconds')

    
    dummy:DRPxPyInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/interactionModelV.pickle.gz')
    dummy.correctFDR(numOfCores=38)
    dummy.pValsHistogram('pValsHistogramCorumDR.png')

    for col in ['interactionPValue', 'interactorPValue', 'XPValue']:
        # dummy.volcanoPlot(f'volcanoPlotModelVCorum{col}DR.png', col, extraFeatures = True)
        dummy.getTopTable(10, col, filepath = f'topTableCorum{col}DR.csv')
    
    dummy.write()

    dummy:DRPxPyInteractionPxModel = read(PATH + '/internal/geneInteractionModel/GLSPValueVAEProteomicsCorum1FDRless0.01/interactionModelV.pickle.gz')
    dummy.correctFDR(numOfCores=38)
    dummy.pValsHistogram('pValsHistogramCorumGeneDependency.png')

    for col in ['interactionPValue', 'interactorPValue', 'XPValue']:
        dummy.volcanoPlot(f'volcanoPlotModelVCorum{col}GeneDependency.png', col, extraFeatures = True)
        dummy.getTopTable(30, col, filepath = f'topTableCorum{col}GeneDependency.csv')
    
    dummy.write()

    dummy:DRPxPyInteractionPxModel = read(PATH + '/internal/interactionModel/String900orBiogrid/interactionModelV.pickle.gz')
    dummy.correctFDR(numOfCores=15)
    dummy.pValsHistogram('pValsHistogramString900andBiogridDR.png')
    
    for col in ['interactionPValue', 'interactorPValue', 'XPValue']:
        # dummy.volcanoPlot(f'volcanoPlotModelVString900andBiogrid{col}DR.png', col, extraFeatures = True)
        dummy.getTopTable(10, col, filepath = f'topTableString900orBiogrid{col}DR.csv')
    
    dummy.write()

    dummy:DRPxPyInteractionPxModel = read(PATH + '/internal/geneInteractionModel/String900orBiogrid/interactionModelV.pickle.gz')
    dummy.correctFDR(numOfCores=11)
    dummy.pValsHistogram('pValsHistogramString900andBiogridGeneDependency.png')
    
    for col in ['interactionPValue', 'interactorPValue', 'XPValue']:
        dummy.volcanoPlot(f'volcanoPlotModelVString900andBiogrid{col}GeneDependecy.png', col, extraFeatures = True)
        dummy.getTopTable(30, col, filepath = f'topTableString900orBiogrid{col}GeneDependency.csv')
    
    dummy.write()




    # dummy.triangulate(-0.06, 0.06, 35, 45, 'Drug Response', 100, 'test.png', True)
    # drugRes.data = drugRes.data.T
    # dummy.scatterTheTopVolcano('topVolcanoPlotScatter.png', ogProteomics, drugRes.data, topNumber=10)

