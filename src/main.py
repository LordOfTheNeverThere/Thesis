# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as t
from statsmodels.stats.multitest import multipletests
from resources import GeneDependency, DRPxPyInteractionPxModel, PyPxDrugInteractionModel, ResidualsLinearModel, ResiduesMatrix, read, PATH, ProteinsMatrix, PairwiseCorrMatrix




if __name__ == '__main__':


    # vaev2 = pd.read_csv(PATH + '/internal/20231023_092657_imputed_proteomics.csv.gz', compression='gzip', index_col=0)
    # #create proteomics matrix
    # proteins = ProteinsMatrix(data = vaev2)
    # #save it 
    # proteins.write(PATH + '/internal/proteomics/proteomicsVAEv2.pickle.gz')
    #create pairwise correlation matrix Using Pearson and GLM
    # proteins:ProteinsMatrix = read(PATH + '/internal/proteomics/proteomicsVAEv2.pickle.gz')
    # pearsonPairwiseVAEv2 = proteins.pearsonCorrelations('coef', 'VAEv2')
    # glmPairwiseVAEv2 = proteins.getGLSCorr("VAEv2")
    # #save them
    # pearsonPairwiseVAEv2.write(PATH + '/internal/pairwiseCorrs/VAEv2/pearsonPairCorr.pickle.gz')
    # glmPairwiseVAEv2.write(PATH + '/internal/pairwiseCorrs/VAEv2/glsPairCorr.pickle.gz')
    # load them
    pearsonPairwiseVAEv2:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAEv2/pearsonPairCorr.pickle.gz')
    glmPairwiseVAEv2:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAEv2/glsPairCorr.pickle.gz')
    #Add all ground truth
    # PairwiseCorrMatrix.addGroundTruths([pearsonPairwiseVAEv2, glmPairwiseVAEv2])
    # pearsonPairwiseVAEv2.aucsCalculator("VAEv2", ["p-value", "coef"], [True, False]) 
    # glmPairwiseVAEv2.aucsCalculator("VAEv2", ["p-value", "coef"], [True, False])
    
    #load all pairwise correlation matrices

    # pearsonPairwiseVAEv2:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAEv2/pearsonPairCorr.pickle.gz')
    # glmPairwiseVAEv2:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAEv2/glsPairCorr.pickle.gz')
    #     #og
    ogPearson:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/OG/baseModelFiltered.pickle.gz')
    #     #Mean
    # # pv75MeanPearson:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr75PV.pickle.gz')
    # # pv80MeanPearson:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr80PV.pickle.gz')
    # # pv75MeanGLM:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz')
    # # pv80MeanGLM:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr80PV.pickle.gz')
    # #     #Vae
    # # vaePearson:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/pearsonPairCorr.pickle.gz')
    # # vaeGLM:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')

    instances = [ogPearson, pearsonPairwiseVAEv2, glmPairwiseVAEv2]
    #Compare all aucs
    aucData = PairwiseCorrMatrix.glsVSPearsonAUC(None, None, "allAUCsBarPlotv2.png", "allAUCs.csv")

    aucData.to_csv('allAUCs.csv', index=False)


