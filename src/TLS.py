import matplotlib.pyplot as plt
from resources import ProteinsMatrix, read, PATH, PairwiseCorrMatrix
import multiprocessing as mp


CPUS = 1
assert CPUS < mp.cpu_count() - 1

if __name__ == '__main__':

    og: ProteinsMatrix = read(PATH + '/datasetsTese/ogProteomics.pickle.gz')
    glsPairwise : PairwiseCorrMatrix = read(PATH + '/datasetsTese/glsPairwiseCorr.pickle.gz')
    glsPairwise = glsPairwise.data
    ppis =  glsPairwise.loc[glsPairwise['corum'] == 1].loc[glsPairwise['beta'] > 0.65].index 



    tlsResidualsMatrix = og.tlsResidues(ppis)
    numOfCells = tlsResidualsMatrix.data.shape[0]  * tlsResidualsMatrix.data.shape[1] # Number of cells in the residuals matrix.

