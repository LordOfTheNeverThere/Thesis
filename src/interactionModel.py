from resources import PATH, read, DRInteractionPxModel, PairwiseCorrMatrix, ProteinsMatrix
import time as t
import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
import matplotlib.pyplot as plt


        

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

    # pca, pcFactors = vaeProteomics.PCA(factorsName='PC', numPC=10)

    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.width', None)
    # pd.set_option('display.max_colwidth', -1)


    dummy = DRInteractionPxModel(ppisOfInterest, ogProteomics, drugRes.data, growthProps, isDrugResSmall = False)
    start = t.time()
    fit = dummy.fit(numOfCores = 38)
    dummy.filepath = PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugLargeRegressor.pickle.gz'
    dummy.write()
    print(f'fitting took {t.time() - start} seconds')

    dummy = DRInteractionPxModel(ppisOfInterest, ogProteomics, drugRes.data, growthProps, isDrugResSmall = True)
    start = t.time()
    fit = dummy.fit(numOfCores = 38)
    dummy.filepath = PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugSmallRegressor.pickle.gz'
    dummy.write()
    print(f'fitting took {t.time() - start} seconds')




    
    # dummy:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugLargeRegressor.pickle.gz')
    # dummy.volcanoPlot('volcanoPlotDrInteractionPxModelDrugLargeExtraSS.png', extraFeatures=False, useExtraSS=False) 
    # dummy.triangulate(-0.06, 0.06, 35, 45, 'Drug Response', 100, 'test.png', True)
    # drugRes.data = drugRes.data.T
    # dummy.scatterTheTopVolcano('topVolcanoPlotScatter.png', ogProteomics, drugRes.data, topNumber=10)



    # # Calculate the effect size of the factor {drug} in linear model on the model's residuals, for small and large model, for all drugs
    # drugSmall:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugSmallRegressor.pickle.gz')
    # start = t.time()
    # drugSmall.resiCorr()
    # print(f'drugSmall.resiCorr() took {t.time() - start} seconds')
    # drugSmall.write()
    # print('done')

    # drugLarge:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugLargeRegressor.pickle.gz')
    # start = t.time()
    # drugLarge.resiCorr()
    # print(f'drugLarge.resiCorr() took {t.time() - start} seconds')
    # drugLarge.write()
    # print('done')

    # drugSmallPCA:DRInteractionPxModel = read(PATH + '/internal/interactionModel/GLPPValueVAEProteomicsCorum1FDRless0.01/drugSmallPCARegressor.pickle.gz')
    # start = t.time()
    # drugSmallPCA.resiCorr()
    # print(f'drugSmallPCA.resiCorr() took {t.time() - start} seconds')