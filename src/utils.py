#%%
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

PATH = "../data"


def getPairwiseCorrData(data: pd.DataFrame, columnName :str ='correlation') -> pd.Series:
    """Gets a series of a pairwise correlations from all pairs of proteins

    Args:
        data (pd.DataFrame): Correlation matrix of proteins x proteins
        columnName (str, optional): The name to give to the column of the series bearing the correlation
         pairwise. Defaults to 'correlation'.

    Returns:
        pd.Series: Series with all pairwise correlations
    """

    data = data.copy()
    row = 0
    col = 1
    columns = data.columns  # list of Columns/Proteins
    proteinIds = []
    pairwiseCorrs = []

    while (col < len(columns)):  # We go every col

        proteinA = columns[row]
        proteinB = columns[col]  # identifier for each interaction
        proteinId = proteinA + ';' + proteinB
        proteinIds.append(proteinId)
        correlation = data.iat[row, col]
        pairwiseCorrs.append(round(correlation,4))

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

# def getGeneIDsCol(data: pd.DataFrame) -> pd.DataFrame:

#     data = data.copy()
#     proteinSubunitsList = [proteinComplex.split(';')
#                            for proteinComplex in list(data.index)]

#     proteinAList = []
#     proteinBList = []

#     for proteinSubunits in proteinSubunitsList:

#         # proteinA = fetchGeneEntrezID(proteinSubunits[0]) # Convert string to Entrez Gene ID
#         # proteinB = fetchGeneEntrezID(proteinSubunits[1])
#         # No conversion need as of now we are working with gene symbols

#         proteinA = proteinSubunits[0]
#         proteinB = proteinSubunits[1]

#         proteinAList.append(proteinA)
#         proteinBList.append(proteinB)

#     # Assigning new cols with the subunit information
#     data['proteinA'] = proteinAList
#     data['proteinB'] = proteinBList

#     return data


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

    pairwiseCorrData.to_csv(PATH + '/datasetsTese/' + fileName + '.csv')

    return pairwiseCorrData


def addGroundTruth(listOfsets: list, data: pd.DataFrame, externalDatasetName: str, filename:str):
    """Append the binary values of a putative PPI, from an external dataset (e.g Corum), to our pairwise correlation Dataframe

    Args:
        listOfsets (list): List of sets of protein protein interactions[{a,b,c},{b,c,d}]
        data (pd.DataFrame): Pairwise correlation dataframe
        externalDatasetName (str): Name to give to the binary column holding the truth value of an PPI is seen in that external Dataset
        filename (str): The name to give to the final file, with the added binary column

    Returns:
        _type_: Data with added column
    """
    data = data.copy()
    data[externalDatasetName] = None

    def addExternalTrueY(model):

        found = 0
        index = 0
        [proteinA, proteinB] = model.name.split(';')
        

        while not found and index < len(listOfsets):

            subset = listOfsets[index]
            if (proteinA in subset) and (proteinB in subset):
                if found:
                    print('Incorrect Code!')
                found = 1
                
                # outputList = outputList.append(1)  # There is interaction

            index += 1
            
        model[externalDatasetName] = found

        return model[externalDatasetName]

    data[externalDatasetName] = data.apply(
        axis=1, func= lambda model: addExternalTrueY(model))
    
    data.to_csv(PATH + '/datasetsTese/' + filename + '.csv')
    
    return data

def getCorumListOfInteractions():
    """DEPRECATED"""
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


    corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')
    listOfSets = corumPPI.apply(axis=1, func=lambda interaction: joinGeneNames(interaction))

    return listOfSets




def getModelsByQuery(datasetToQuery: str, featureToQuery:str, valueToQuery)->set:
    """This func does a query on 3 of possible datasets (samplesheet, drugresponse, CRISPR) and returns a list of the models that
    abide by this query

    Args:
        datasetToQuery (str): name of the dataset to query (samplesheet, drugresponse, CRISPR)
        featureToQuery (str): the column name to query in the dataset
        valueToQuery (_type_):the value to query in the feature/column of the dataset

    Returns:
        set: Returns a set composed of all models abiding by the query
    """

    if datasetToQuery == 'samplesheet': # In the case we are querying the samplesheet.csv for model specific features
        datasetToQuery: pd.DataFrame = pd.read_csv(PATH + "/datasetsTese/" + datasetToQuery + ".csv", index_col='model_id')
    elif datasetToQuery == 'drugresponse':
        datasetToQuery: pd.DataFrame = pd.read_csv(
            PATH + "/datasetsTese/" + datasetToQuery + ".csv", index_col='model_id')
    elif datasetToQuery == 'CRISPR':
        datasetToQuery: pd.DataFrame = pd.read_csv(
            PATH + "/datasetsTese/crisprcas9_22Q2.csv", index_col='model_id')
    else:
        try:
            datasetToQuery: pd.DataFrame = pd.read_csv(
                PATH + "/datasetsTese/" + datasetToQuery + ".csv", index_col='model_id')
        except:
            print("The only values accepted for datasetToQuery are: samplesheet, drugresponse, CRISPR. You inserted \n" + datasetToQuery)
        

    queriedDataset = datasetToQuery.query(f"{featureToQuery} == @valueToQuery")

    return set(queriedDataset.index)

def getUniqueSetValues(filepath: str, feature: str):
    """Returns a set of unique values from a feature of a dataframe

    Args:
        filepath (str): Filepath to load the dataframe
        feature (str): The column name to extract the unique set of values

    Returns:
        set: The uniqye set of values in a column of a Dataframe
        dict: Dictionary with the keys are the unique values in a column and the values(of the dict) as the number of
        occurances in of each value(of the feature)
    """

    
    data = pd.read_csv(index_col='model_id', filepath_or_buffer=filepath)
    setOfValues = data[feature].unique()
    setOfValues = set(setOfValues)
    occurancesDict = data.groupby(
        feature).count().to_dict()['model_name']


    return setOfValues, occurancesDict
# %%
