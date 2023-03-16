import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt


def getPairwiseCorrData(data: pd.DataFrame, columnName='correlation') -> pd.Series:

    data = data.copy()
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

    pairwiseCorrsSeries = pd.Series(
        data=pairwiseCorrs, index=proteinIds, name=columnName)
    return pairwiseCorrsSeries


# def fetchGeneEntrezID(proteinName: str, geneIDsData: pd.DataFrame = pd.read_table('../externalDatasets/geneIDsNCBI.tsv')) -> int:

#     geneID = geneIDsData.where(geneIDsData['Symbol'] == proteinName)[
#         'NCBI GeneID']
#     # Clean all the nan values and get the first value which is the desired ID
#     geneID = geneID.dropna()[0]

#     return geneID

# No conversion needed we are working with gene Symbols

def getGeneIDsCol(data: pd.DataFrame) -> pd.DataFrame:

    data = data.copy()
    proteinSubunitsList = [proteinComplex.split('-+-')
                           for proteinComplex in list(data.index)]

    proteinAList = []
    proteinBList = []

    for proteinSubunits in proteinSubunitsList:

        # proteinA = fetchGeneEntrezID(proteinSubunits[0]) # Convert string to Entrez Gene ID
        # proteinB = fetchGeneEntrezID(proteinSubunits[1])
        # No conversion need as of now we are working with gene symbols

        proteinA = proteinSubunits[0]
        proteinB = proteinSubunits[1]

        proteinAList.append(proteinA)
        proteinBList.append(proteinB)

    # Assigning new cols with the subunit information
    data['proteinA'] = proteinAList
    data['proteinB'] = proteinBList

    return data


def getPairwiseCorrelation(data: pd.DataFrame, fileName: str, columnName: str) -> tuple:
    """_summary_

    Args:
        data (pd.DataFrame): Proteomics Data in transposed form
        fileName (str): Name of the csv file created with the pairwise correlation

    Returns:
        tuple: Tuple with Dataframes with the pairwise correlation, the first one with subuints of the protein protein complex
    """
    data = data.copy()

    pearsonCorrMatrix = data.corr(method='pearson')
    pairwiseCorrData = getPairwiseCorrData(
        data=pearsonCorrMatrix, columnName=columnName)
    pairwiseCorrData.sort_values(ascending=False, inplace=True)
    proteinIDCorrData = getGeneIDsCol(data=pairwiseCorrData)

    proteinIDCorrData.to_csv(
        '../data/datasetsTese/' + fileName + '.csv')

    return proteinIDCorrData, pairwiseCorrData


def addGroundTruth(listOfsets: list, data: pd.DataFrame, externalDatasetName: str):
    
    data = data.copy()
    data[externalDatasetName] = None

    def addExternalTrueY(model):

        found = 0
        index = 0

        while not found and index < len(listOfsets):

            subset = listOfsets[index]
            if (model['proteinA'] in subset) and (model['proteinB'] in subset):
                found = 1
                
                # outputList = outputList.append(1)  # There is interaction

            index += 1
            
        model[externalDatasetName] = found

        return model[externalDatasetName]

    data[externalDatasetName] = data.apply(
        axis=1, func= lambda model: addExternalTrueY(model))
    
    return data

def getCorumListOfInteractions():
    def joinGeneNames(interaction):
        subset1 = interaction['subunits(Gene name)']
        subset2 = interaction['subunits(Gene name syn)']
        subset1 = subset1.split(';')
        subset2 = subset2.split(';')
        lenghtSubset1 = len(subset1)

        for index in range(0, lenghtSubset1):
            proteinAliases = subset2[index]
            proteinAliases = proteinAliases.split(' ')
            subset1 = subset1 + proteinAliases

        subset1 = set(subset1)
        subset1.discard('None')
        subset1.discard('')

        return subset1


    corumPPI = pd.read_json('../data/externalDatasets/corumPPI.json')
    listOfSets = corumPPI.apply(axis=1, func=lambda interaction: joinGeneNames(interaction))

    return listOfSets
