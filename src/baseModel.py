# Imports
import pandas as pd
from classes import ProteinsMatrix, PairwiseCorrMatrix, ppiDataset
import pickle
import gzip

from env import PATH






if __name__ == '__main__':

    proteins = ProteinsMatrix(
            PATH + '/datasetsTese/proteomicsDataTrans.csv', index_col='modelID')
    
    

    pairwiseCorrPValues = proteins.pearsonCorrelations('globalCorrelation')
    # External ppi data
    corum = ppiDataset(PATH + '/externalDatasets/corumPPI.csv.gz', compression='gzip')
    biogrid = ppiDataset(PATH + '/externalDatasets/biogridPPIHuman.csv.gz', compression='gzip')
    string = ppiDataset(PATH + '/externalDatasets/stringPPI900Selected.csv.gz', compression='gzip')

    #Add ground truths to the pairwise coor matrix
    for dataset in [(corum,'corum'), (biogrid,'biogrid'), (string,'string')]:
        pairwiseCorrPValues.addGroundTruth(dataset[0].getPPIs(dataset[1]), dataset[1], None)

    pairwiseCorrPValues.aucCalculator('corum', 'baseModel')    

    

    with gzip.open(PATH + '/datasetsTese/'+ 'baseModelPairwiseCorr.pickle.gz', 'wb') as f:
        pickle.dump(pairwiseCorrPValues, f)
    f.close()
    print(pairwiseCorrPValues.data)
    print(pairwiseCorrPValues.auc)


    # Converting .csv.gz to pickle.gz

    #BaseModel

    # baseModel = PairwiseCorrMatrix(PATH+'/datasetsTese/baseModel.csv.gz', data=None, compression='gzip', index_col='PPI')
    # baseModel.data.rename(columns={'Corum': 'corum'}, inplace = True)
    # baseModel.aucCalculator('corum', 'baseModel')
    # print(baseModel.data)
    # print(baseModel.auc)

    # with gzip.open(PATH + '/datasetsTese/'+ 'baseModelPairwiseCorr.pickle.gz', 'wb') as f:
    #     pickle.dump(baseModel, f)

    #MOFA

    # mofa = PairwiseCorrMatrix(PATH + '/datasetsTese/baseMOFAPairwiseCorr.csv.gz', data=None, compression='gzip', index_col='PPI')
    # mofa.aucCalculator('corum', 'mofa')
    # with gzip.open(PATH + '/datasetsTese/'+ 'mofaPairwiseCorr.pickle.gz', 'wb') as f:
    #     pickle.dump(mofa, f)

    # #gls

    # with gzip.open(PATH + '/datasetsTese/' + 'glsPairwiseCorr.pickle.gz', 'rb') as f:
    #     gls = pickle.load(f)
    # f.close()
    # gls.data.sort_values(by='glsCoefficient', ascending=False, inplace=True)
    # gls.aucCalculator('corum', 'gls')
    # print(gls.data)
    # print(gls.auc)
    # with gzip.open(PATH + '/datasetsTese/'+ 'glsPairwiseCorr.pickle.gz', 'wb') as f:
    #     pickle.dump(gls, f)
    # f.close()

    # #Corum

    # corum = ppiDataset(filepath=PATH + '/externalDatasets/corumPPI.csv.gz', compression='gzip')
    # corum.getPPIs('corum')
    # print(corum.ppis)
    # with gzip.open(PATH + '/externalDatasets/'+ 'corum.pickle.gz', 'wb') as f:
    #     pickle.dump(corum, f)

    # #Biogrid

    # biogrid = ppiDataset(filepath=PATH + '/externalDatasets/biogridPPIHuman.csv.gz', compression='gzip')
    # biogrid.getPPIs('biogrid')
    # print(biogrid.ppis)
    # with gzip.open(PATH + '/externalDatasets/'+ 'biogrid.pickle.gz', 'wb') as f:
    #     pickle.dump(biogrid, f)

    # #String

    # string = ppiDataset(filepath=PATH + '/externalDatasets/stringPPI900Selected.csv.gz', compression='gzip')
    # string.getPPIs('string')
    # print(string.ppis)
    # with gzip.open(PATH + '/externalDatasets/'+ 'string.pickle.gz', 'wb') as f:
    #     pickle.dump(string, f)

    # #ogProteomics

    # ogProteomics = ProteinsMatrix(PATH + '/datasetsTese/proteomicsDataTrans.csv')
    # print(ogProteomics.data)
    # with gzip.open(PATH + '/datasetsTese/'+ 'ogProteomics.pickle.gz', 'wb') as f:
    #     pickle.dump(ogProteomics, f)
    
    # # #mofaProteomics

    # mofaProteomics = ProteinsMatrix(PATH + '/datasetsTese/proteomicsMOFA.csv.gz', index_col='Unnamed: 0', compression='gzip')
    # print(mofaProteomics.data)
    # with gzip.open(PATH + '/datasetsTese/' + 'mofaProteomics.pickle.gz', 'wb') as f:
    #     pickle.dump(mofaProteomics, f)

        
    
