import matplotlib.pyplot as plt
from resources import ProteinsMatrix, read, PATH
import multiprocessing as mp


CPUS = 5
assert CPUS < mp.cpu_count() - 1

def wrapper(proteomics:ProteinsMatrix, filepath:str):
    
    PairwiseCorr = proteomics.pearsonCorrelations('pearsonR')
    PairwiseCorr.write(filepath)

og: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
vae= ProteinsMatrix(PATH + '/datasetsTese/proteomicsVAE.csv.gz', compression='gzip', index_col='Unnamed: 0')
vae.data = vae.data.iloc[0:949,:]

filepaths = [PATH + '/datasetsTese/baseModelFiltered.pickle.gz', PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz']
proteomics = [og, vae]



with mp.Pool(CPUS) as process:
    checkPPIGen = process.starmap(wrapper, zip(proteomics, filepaths))  # While Cycle


# glsCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
# pearsonCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/baseModelFiltered.pickle.gz') # R
# vaeGLSCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz') # gls + vae
# vaePearsonCorrs:PairwiseCorrMatrix = read(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz') # R + vae

# print('Loading Completed')

##Add all aucs
# pairwiseCorrs = [glsCorrs, pearsonCorrs, vaeGLSCorrs, vaePearsonCorrs]
# yColumnLists = [['corum', 'corum'],['corum', 'corum'],['corum', 'corum'],['corum', 'corum']]
# proxyColumnLists = [['pValue', 'beta'],['pValue', 'pearsonR'],['pValue', 'beta'],['pValue', 'pearsonR']]
# ascendingLists = [[True, False],[True, 'pearsonR'],[True, False],[True, False]]
# labels = ['gls', 'pearson', 'VAE-GLS', 'VAE-pearson']
# filepaths= [PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz', PATH + '/datasetsTese/baseModelFiltered.pickle.gz', PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz', PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz']




# glsCorrs.aucCalculator('corum', 'gls','pValue', True)
# glsCorrs.aucCalculator('corum', 'gls','beta', False)
# print(glsCorrs)
# glsCorrs.write(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')

# pearsonCorrs.aucCalculator('corum', 'pearson','pValue', True)
# pearsonCorrs.aucCalculator('corum', 'pearson','beta', False)
# print(pearsonCorrs)
# pearsonCorrs.write(PATH + '/datasetsTese/baseModelFiltered.pickle.gz')

# vaeGLSCorrs.aucCalculator('corum', 'VAE-GLS','pValue', True)
# vaeGLSCorrs.aucCalculator('corum', 'VAE-GLS','beta', False)
# print(vaeGLSCorrs)
# vaeGLSCorrs.write(PATH + '/datasetsTese/VAEGLSPairCorr.pickle.gz')

# vaePearsonCorrs.aucCalculator('corum', 'VAE-pearson','pValue', True)
# vaePearsonCorrs.aucCalculator('corum', 'VAE-pearson','beta', False)
# print(vaePearsonCorrs)
# vaePearsonCorrs.write(PATH + '/datasetsTese/VAEPearsonPairCorr.pickle.gz')





# tlsResidualsMatrix = proteomics.tlsResidues(glsPairwise)
# print(glsPairwise.data.loc[glsPairwise.data['p-value'] > 0])
# print(tlsResidualsMatrix.data.shape)
# print(tlsResidualsMatrix.data.isna().sum().sum())

