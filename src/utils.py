import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt


def getPairwiseCorrData(data: pd.DataFrame) -> pd.Series:

    row = 0
    col = 1
    columns = data.columns  # list of Columns/Proteins
    proteinIds = []
    pairwiseCorrs = []

    while (col < len(columns)):  # We go every col

        proteinA = columns[row]
        proteinB = columns[col]  # identifier for each interaction
        proteinId = proteinA + '-+-' + proteinB
        proteinIds.append(proteinId)
        correlation = data.iat[row, col]
        pairwiseCorrs.append(correlation)

        row += 1
        if row == col:  # And then every row, only going to the next col when the number of rows catches up to the col number, so we would be in the diagonal of the corr matric
            col += 1
            row = 0

    pairwiseCorrsSeries = pd.Series(data=pairwiseCorrs, index=proteinIds)
    return pairwiseCorrsSeries


def fetchGeneEntrezID(proteinName: str, geneIDsData: pd.DataFrame = pd.read_table('../externalDatasets/geneIDsNCBI.tsv')) -> int:

    geneID = geneIDsData.where(geneIDsData['Symbol'] == proteinName)[
        'NCBI GeneID']
    # Clean all the nan values and get the first value which is the desired ID
    geneID = geneID.dropna()[0]

    return geneID


def getGeneIDsCol(data: pd.DataFrame) -> pd.DataFrame:

    proteinSubunitsList = [proteinComplex.split('-+-')
                           for proteinComplex in list(data.index)]
    
    proteinAList = proteinBList = proteinABList = []
    
    for proteinSubunits in proteinSubunitsList:

        proteinA = fetchGeneEntrezID(proteinSubunits[0]) # Convert string to Entrez Gene ID
        proteinB = fetchGeneEntrezID(proteinSubunits[1])

        proteinAList.append(proteinA)
        proteinBList.append(proteinB)
        proteinABList.append(proteinA + ',' + proteinB)
    
    #Assigning new cols with the subunit information
    data['proteinA'] = proteinAList; data['proteinB'] = proteinBList; data['proteinAB'] = proteinABList

    return data 


    
