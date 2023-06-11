# Imports
import pandas as pd
import matplotlib.pyplot as plt
import time as t
from resources import ResidualsLinearModel, read, PATH








if __name__ == '__main__':

    drugRes = read(PATH + '/datasetsTese/drugResponse.pickle.gz')

    Y: pd.DataFrame = drugRes.data.T # Samples should be rows and not columns
    Y = Y.fillna(Y.mean(axis=0)) 

    residuals = read(PATH + '/datasetsTese/residualsMatrixGLSGreater0.65.pickle.gz')

    X = residuals.data

    samplesheet = pd.read_csv(PATH + '/datasetsTese/samplesheet.csv', index_col='model_id')
    otherCovariates = samplesheet[['tissue', 'growth_properties']].dropna(axis=0, how='any')
    otherCovariates['hymCellLine'] = (otherCovariates['tissue'] == 'Haematopoietic and Lymphoid').astype(int)
    otherCovariates['lung'] = (otherCovariates['tissue'] == 'lung').astype(int)
    otherCovariates = pd.get_dummies(otherCovariates, columns=['growth_properties'], prefix='', prefix_sep='')
    otherCovariates = otherCovariates.drop(columns=['tissue'])


    regressor = ResidualsLinearModel(Y, X, otherCovariates, residualsType='malahanobis')
    print(regressor.fit_matrix())
    regressor.volcanoPlot(filepath='testVolcano.png')
    