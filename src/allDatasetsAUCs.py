import pandas as pd
import seaborn as sb
import numpy as np
from classes import ppiDataset
from env import PATH





if __name__ == "__main__":

    biogrid = pd.read_table(PATH + '/externalDatasets/biogridPPIALL.txt')
    biogrid = biogrid.query("`Organism Name Interactor A` == 'Homo sapiens' & `Organism Name Interactor B` == 'Homo sapiens'")
    biogrid.to_csv(PATH + '/externalDatasets/biogridPPIHuman.csv.gz', compression='gzip')
    
