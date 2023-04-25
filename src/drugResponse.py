import pandas as pd
from classes import DrugResponseMatrix
from env import PATH




if __name__ == "__main__":

    drugResponse = pd.read_csv(PATH + '/datasetsTese/drugResponse.csv', index_col='Unnamed: 0')
    maxConcentration = pd.read_csv(PATH + '/externalDatasets/drugMaxScreenConcentration.csv', index_col='Unnamed: 0')


    print(tuple([maxConcentration['MAX_CONC'],maxConcentration.index]))