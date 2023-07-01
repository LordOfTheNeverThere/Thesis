import matplotlib.pyplot as plt
from resources import ProteinsMatrix, read, PATH, PairwiseCorrMatrix
import multiprocessing as mp


if __name__ == '__main__':

    #Load all pairwise correlation matrices

    # og:PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/pearsonPairCorr.pickle.gz')
    # og.filepath = PATH + '/internal/pairwiseCorrs/VAE/pearsonPairCorr.pickle.gz'
    # og.corrCumSums = {}
    # og.indexes = {}
    # og.aucs = {}
    # og.labels = {}
    # og.yColumn = 'corum'
    # og.proxies = ['pearsonR', 'pValue']
    # og.ascendings = [False, True]
    # og.proteomicsType = 'VAE'
    # PairwiseCorrMatrix.getAucs([og])




    proteinMean75PV:ProteinsMatrix = read(PATH + '/internal/proteomics/mean75PVProteomics.pickle.gz')
    proteinMean80PV:ProteinsMatrix = read(PATH + '/internal/proteomics/mean80PVProteomics.pickle.gz')

    proteinMean75Pair: PairwiseCorrMatrix = proteinMean75PV.pearsonCorrelations('pearsonR', '75PVMean')
    proteinMean80Pair: PairwiseCorrMatrix = proteinMean80PV.pearsonCorrelations('pearsonR', '80PVMean')

    proteinMean75Pair.filepath = PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr75PV.pickle.gz'
    proteinMean80Pair.filepath = PATH + '/internal/pairwiseCorrs/Mean/pearsonPairCorr80PV.pickle.gz'

    proteinMean75Pair.write()
    proteinMean80Pair.write()

    proteinMean80Pair = proteinMean80PV.getGLSCorr('80PVMean', coefColumnName='beta')
    proteinMean75Pair = proteinMean75PV.getGLSCorr('75PVMean', coefColumnName='beta')

    proteinMean75Pair.filepath = PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr75PV.pickle.gz'
    proteinMean80Pair.filepath = PATH + '/internal/pairwiseCorrs/Mean/glsPairCorr80PV.pickle.gz'

    proteinMean75Pair.write()
    proteinMean80Pair.write()