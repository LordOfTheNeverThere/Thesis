#%%
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from classes import PairwiseCorrMatrix

from env import PATH


# def getPairwiseCorrData(data: pd.DataFrame, columnName :str ='correlation') -> pd.Series:
#     """ DEPRECATED
#     Gets a series of a pairwise correlations from all pairs of proteins

#     Args:
#         data (pd.DataFrame): Correlation matrix of proteins x proteins
#         columnName (str, optional): The name to give to the column of the series bearing the correlation
#          pairwise. Defaults to 'correlation'.

#     Returns:
#         pd.Series: Series with all pairwise correlations
#     """

#     data = data.copy()
#     row = 0
#     col = 1
#     columns = data.columns  # list of Columns/Proteins
#     proteinIds = []
#     pairwiseCorrs = []

#     while (col < len(columns)):  # We go every col

#         proteinA = columns[row]
#         proteinB = columns[col]  # identifier for each interaction
#         proteinId = proteinA + ';' + proteinB
#         proteinIds.append(proteinId)
#         correlation = data.iat[row, col]
#         pairwiseCorrs.append(round(correlation,4))

#         row += 1
#         if row == col:  # And then every row, only going to the next col when the number of rows catches up to the col number, so we would be in the diagonal of the corr matric
#             col += 1
#             row = 0

#     pairwiseCorrsSeries = pd.Series(
#         data=pairwiseCorrs, index=proteinIds, name=columnName)
#     return pairwiseCorrsSeries


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


def getPairwiseCorrelation(data: pd.DataFrame, fileName: str, columnName: str, counting:bool = True) -> pd.DataFrame:
    """_summary_

    Args:
        data (pd.DataFrame): Proteomics Data in transposed form
        fileName (str): Name of the csv file created with the pairwise correlation

    Returns:
        Dataframe: Dataframe with the pairwise correlation
    """
    data = data.copy()
    # Correlation Matrix
    pearsonCorrMatrix = data.corr(method='pearson').round(decimals=4)
    
    pairwiseCorrData = pearsonCorrMatrix.where(np.triu(np.ones(pearsonCorrMatrix.shape), k=1).astype(bool)).stack().reset_index()
    pairwiseCorrData['PPI'] = pairwiseCorrData['level_0'] +";" + pairwiseCorrData['level_1']
    pairwiseCorrData.drop(columns=['level_0', 'level_1'], inplace=True)
    pairwiseCorrData = pairwiseCorrData.set_index('PPI')
    pairwiseCorrData.columns=[columnName]
    
    if counting:
        
        # Co-occorance Matrix
        coOccuranceMatrix = (data/data).fillna(0).astype(int) #Get 1 and 0 depending on if there is a value or not in a specific spot
        coOccuranceMatrix = coOccuranceMatrix.T.dot(coOccuranceMatrix)  #Simple linear algebra to get the co-occurance values
        coOccuranceData = coOccuranceMatrix.where(np.triu( np.ones(coOccuranceMatrix.shape), k=1).astype(bool)).stack().reset_index()
        coOccuranceData['PPI'] = coOccuranceData['level_0'] +";" + coOccuranceData['level_1']
        coOccuranceData.drop(columns=['level_0', 'level_1'], inplace=True)
        coOccuranceData = coOccuranceData.set_index('PPI')
        coOccuranceData.columns = ['counts']
        pairwiseCorrData = pairwiseCorrData.merge(
            coOccuranceData, on='PPI', how='left')


    
    
    pairwiseCorrData.sort_values(by=columnName,ascending=False, inplace=True)
    
    if fileName:
        pairwiseCorrData.to_csv(PATH + '/datasetsTese/' + fileName + '.csv')

    return pairwiseCorrData


def addGroundTruth(ppis: set, data: pd.DataFrame, externalDatasetName: str, filename:str = None):
    """Append the binary values of a putative PPI, from an external dataset (e.g Corum), to our pairwise correlation Dataframe

    Args:
        ppis (set): sets of tuples of protein protein interactions{(a,b),(c,b),(c,d)}
        data (pd.DataFrame): Pairwise correlation dataframe
        externalDatasetName (str): Name to give to the binary column holding the truth value of an PPI is seen in that external Dataset
        filename (str): The name to give to the final file, with the added binary column

    Returns:
        _type_: Data with added column
    """
    data = data.copy()

    def addExternalTrueY(model):

        found = 0
        [proteinA, proteinB] = model.name.split(';')
        ppiAB: tuple = (proteinA, proteinB) #In my implementation the ppis have (A,B) but not (B,A), they are combinations
        ppiBA: tuple = (proteinB, proteinA)
        
        if ppiAB in ppis or ppiBA in ppis:
            found = 1
            
        model[externalDatasetName] = found

        return model[externalDatasetName]

    data[externalDatasetName] = data.apply(
        axis=1, func= lambda model: addExternalTrueY(model))
    
    if filename:

        data.to_csv(PATH + '/datasetsTese/' + filename + '.csv')
    
    return data


# DEPRECATED
# def addGroundTruthTreeNode(ppiTree: TreeNode, data: pd.DataFrame, externalDatasetName: str, filename:str = None):
#     """Append the binary values of a putative PPI, from an external dataset (e.g Corum), to our pairwise correlation Dataframe

#     Args:
#         listOfsets (list): List of sets of protein protein interactions[{a,b,c},{b,c,d}]
#         data (pd.DataFrame): Pairwise correlation dataframe
#         externalDatasetName (str): Name to give to the binary column holding the truth value of an PPI is seen in that external Dataset
#         filename (str): The name to give to the final file, with the added binary column

#     Returns:
#         _type_: Data with added column
#     """
#     data = data.copy()
#     data[externalDatasetName] = None

#     def addExternalTrueY(model):

#         found = 0
#         index = 0
#         [proteinA, proteinB] = model.name.split(';')
        

#         proteinANode :TreeNode  = ppiTree.getNodeFirstLayer(proteinA)
        
#         if proteinANode and proteinB in proteinANode.getChildrenValue():
#             found = 1
            
#         model[externalDatasetName] = found

#         return model[externalDatasetName]

#     data[externalDatasetName] = data.apply(
#         axis=1, func= lambda model: addExternalTrueY(model))
    
#     if filename:

#         data.to_csv(PATH + '/datasetsTese/' + filename + '.csv')
    
#     return data


# def getCorumListOfInteractions():
#     """DEPRECATED"""
#     def joinGeneNames(interaction):
#         subset1 = interaction['subunits(Gene name)']
#         subset2 = interaction['subunits(Gene name syn)']
#         subset1 = subset1.split(';')
#         subset2 = subset2.split(';')
#         lenghtSubset1 = len(subset1)

#         for index in range(0, lenghtSubset1):
#             proteinAliases = subset2[index]
#             proteinAliases = proteinAliases.split(' ')
#             subset1 = subset1 + proteinAliases

#         subset1 = set(subset1)
#         subset1.discard('None')
#         subset1.discard('')

#         return subset1


#     corumPPI = pd.read_json(PATH + '/externalDatasets/corumPPI.json')
#     listOfSets = corumPPI.apply(axis=1, func=lambda interaction: joinGeneNames(interaction))

#     return listOfSets




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

def drawRecallCurves(paiwiseMatrices : list[PairwiseCorrMatrix], colours: list, filename: str):

    _, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)

    for index, pairwiseCorr in enumerate(paiwiseMatrices):
        ax.plot(
            pairwiseCorr.indexes,
            pairwiseCorr.corrCumSum,
            label=pairwiseCorr.label,
            c=colours[index],
        )

    ax.plot([0, 1], [0, 1], "k--", lw=0.3)
    ax.legend(loc="lower right", frameon=False)

    ax.set_ylabel("Cumulative sum")
    ax.set_xlabel("Ranked correlation")
    ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

    plt.savefig(PATH +  '/images/' + filename, bbox_inches="tight")
    plt.close("all")

