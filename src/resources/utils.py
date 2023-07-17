
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip
from resources import PATH, MatrixData, PairwiseCorrMatrix


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


# def fetchGeneEntrezID(proteinName: str, geneIDsData: pd.DataFrame = pd.read_table('../external/geneIDsNCBI.tsv')) -> int:

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


def getPairwiseCorrelation(data: pd.DataFrame, fileName: str, columnName: str, counting: bool = True) -> pd.DataFrame:
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

    pairwiseCorrData = pearsonCorrMatrix.where(np.triu(
        np.ones(pearsonCorrMatrix.shape), k=1).astype(bool)).stack().reset_index()
    pairwiseCorrData['PPI'] = pairwiseCorrData['level_0'] + \
        ";" + pairwiseCorrData['level_1']
    pairwiseCorrData.drop(columns=['level_0', 'level_1'], inplace=True)
    pairwiseCorrData = pairwiseCorrData.set_index('PPI')
    pairwiseCorrData.columns = [columnName]

    if counting:

        # Co-occorance Matrix
        # Get 1 and 0 depending on if there is a value or not in a specific spot
        coOccuranceMatrix = (data/data).fillna(0).astype(int)
        # Simple linear algebra to get the co-occurance values
        coOccuranceMatrix = coOccuranceMatrix.T.dot(coOccuranceMatrix)
        coOccuranceData = coOccuranceMatrix.where(np.triu(
            np.ones(coOccuranceMatrix.shape), k=1).astype(bool)).stack().reset_index()
        coOccuranceData['PPI'] = coOccuranceData['level_0'] + \
            ";" + coOccuranceData['level_1']
        coOccuranceData.drop(columns=['level_0', 'level_1'], inplace=True)
        coOccuranceData = coOccuranceData.set_index('PPI')
        coOccuranceData.columns = ['counts']
        pairwiseCorrData = pairwiseCorrData.merge(
            coOccuranceData, on='PPI', how='left')

    pairwiseCorrData.sort_values(by=columnName, ascending=False, inplace=True)

    if fileName:
        pairwiseCorrData.to_csv(PATH + '/internal/' + fileName + '.csv')

    return pairwiseCorrData


def addGroundTruth(ppis: set, data: pd.DataFrame, externalDatasetName: str, filename: str = None):
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
        # In my implementation the ppis have (A,B) but not (B,A), they are combinations
        ppiAB: tuple = (proteinA, proteinB)
        ppiBA: tuple = (proteinB, proteinA)

        if ppiAB in ppis or ppiBA in ppis:
            found = 1

        model[externalDatasetName] = found

        return model[externalDatasetName]

    data[externalDatasetName] = data.apply(
        axis=1, func=lambda model: addExternalTrueY(model))

    if filename:

        data.to_csv(PATH + '/internal/' + filename + '.csv')

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

#         data.to_csv(PATH + '/internal/' + filename + '.csv')

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


#     corumPPI = pd.read_json(PATH + '/external/corumPPI.json')
#     listOfSets = corumPPI.apply(axis=1, func=lambda interaction: joinGeneNames(interaction))

#     return listOfSets


def getModelsByQuery(datasetToQuery: MatrixData, featureToQuery: str, valueToQuery) -> set:
    """This func does a query on 3 of possible datasets (samplesheet, drugresponse, CRISPR) and returns a list of the models that
    abide by this query

    Args:
        datasetToQuery (str): name of the dataset to query (samplesheet, drugresponse, CRISPR)
        featureToQuery (str): the column name to query in the dataset
        valueToQuery (_type_):the value to query in the feature/column of the dataset

    Returns:
        set: Returns a set composed of all models abiding by the query
    """

    queriedDataset = datasetToQuery.query(f"{featureToQuery} == @valueToQuery")

    return set(queriedDataset.index)


def getUniqueSetValues(data: pd.DataFrame, feature: str):
    """Returns a set of unique values from a feature of a dataframe

    Args:
        filepath (str): Filepath to load the dataframe
        feature (str): The column name to extract the unique set of values

    Returns:
        set: The uniqye set of values in a column of a Dataframe
        dict: Dictionary with the keys are the unique values in a column and the values(of the dict) as the number of
        occurances in of each value(of the feature)
    """

    setOfValues = data[feature].unique()
    setOfValues = set(setOfValues)
    occurancesDict = data.groupby(
        feature).count().to_dict()[feature]

    return setOfValues, occurancesDict


def drawRecallCurves(paiwiseMatrices: list[PairwiseCorrMatrix], colours: list, filename: str, proxyColumn: str, yColumn: str):

    _, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600)

    for index, pairwiseCorr in enumerate(paiwiseMatrices):
        ax.plot(
            pairwiseCorr.indexes[proxyColumn][yColumn],
            pairwiseCorr.corrCumSums[proxyColumn][yColumn],
            label=pairwiseCorr.labels[proxyColumn][yColumn],
            c=colours[index],
        )

    ax.plot([0, 1], [0, 1], "k--", lw=0.3)
    ax.legend(loc="lower right", frameon=False)

    ax.set_ylabel("Cumulative sum")
    ax.set_xlabel("Ranked correlation")
    ax.grid(True, axis="both", ls="-", lw=0.1, alpha=1.0, zorder=0)

    plt.savefig('../images/' + filename, bbox_inches="tight")
    plt.close("all")


def read(filepath: str):
    """Load one of the pickled objects stored in filepath

    Args:
        filepath (str): filepath of pickled, gziped object

    Returns:
        _type_: object
    """
    import sys
    sys.path.append('resources')
    with gzip.open(filepath, 'rb') as f:
        object = pickle.load(f)
    f.close()

    object.filepath = filepath

    return object


def pxPyScatterPlots(other: PairwiseCorrMatrix, limitsR: tuple[float], limitsMetricOther: tuple[ﬂoat], corum: int, lung: bool, sortingColum: str, ascending: bool) -> None:
    """Creates a scatter plot with top 5 PPI according with queries done to the 

    Args:
        limitsMetric (tuple[float]): Interval (start, end) for metric for
        limitsMetricOther (tuple[): _description_
        corum (int): _description_
        lung (bool): _description_
        sortingColum (str): _description_
        ascending (bool): _description_
    """

    # Dissect Code And Restructure

    """    limitsR = (-0.01,0.01)
    limitsGLS = (-0.5,0.5)
    corum = 1
    lung = True
    sortingColumn = 'glsCoefficient'
    ascending = False
    sortingSymbol = '↑' if ascending else '↓'


    baseModel :PairwiseCorrMatrix = read(PATH + '/internal/baseModelFiltered.pickle.gz')
    glsPairwise :PairwiseCorrMatrix = read(PATH + '/internal/glsPairwiseCorr.pickle.gz')
    ogProteomics :ProteinsMatrix = read(PATH + '/internal/ogProteomics.pickle.gz')
    sampleSheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col='model_id')

    # Merging/sorting and quering

    mergedData:pd.DataFrame = baseModel.data.merge(right = glsPairwise.data[['glsCoefficient']], on ='PPI')
    mergedData = mergedData.loc[mergedData['corum'] == corum]
    # mergedData = mergedData.loc[mergedData['biogrid'] == 0]
    # mergedData = mergedData.loc[mergedData['string'] == 0]
    # mergedData = mergedData.loc[mergedData['string400'] == 0]

    
    mergedData = mergedData.query('globalCorrelation < @limitsR[1] & globalCorrelation > @limitsR[0] & glsCoefficient > @limitsGLS[0] & glsCoefficient < @limitsGLS[1]').copy().sort_values([sortingColumn], ascending=[ascending])
    print(mergedData)
    
    top5evidencePPI = mergedData.head(5)
    indexes = list(top5evidencePPI.index)
    setOfPPIs = [(proteins.split(';')[0], proteins.split(';')[1] ) for proteins in indexes]#Unpack PPIs of opposing inensities into tuples of proteins 
    proteins =  [protein for sublist in setOfPPIs for protein in sublist]
    print(setOfPPIs, proteins)
    fig, ax = plt.subplots(10, 1, figsize=(20, 70))
    for index, protein in enumerate(proteins):
        
        proteinExpression =  ogProteomics.data.loc[:,protein] #Get the protein expresion of the two proteins in the ppi
        proteinExpression = proteinExpression.dropna(axis=0)

        samples = list(proteinExpression.index) #Get all samples where both proteins have an expression value
        mVSamples = [ogProteomics.data.loc[sample,:].isna().sum()  for sample in samples] #Count the #MV in each samples
        meanValsSamples = [ogProteomics.data.loc[sample,:].mean()  for sample in samples] #Count the #MV in each samples
        # replicateCountsSamples = [sampleSheet.loc[sample,'replicate_correlation']  for sample in samples] #Count the #MV in each samples
        # print(replicateCountsSamples)

        if not lung: # If we dont want a specific tisssue like lung

            samples = [sample for sample in samples if sampleSheet.loc[sample, 'tissue']!='Lung']
            proteinAB = proteinAB.loc[samples,:]
        
        tissues = [sampleSheet.loc[sample, 'tissue'] for sample in samples] #get all sample's tissues
        



        tissues = [TISSUEPALETTE[tissue] for tissue in tissues] #convert tissues to their respective colours


        ax[index].scatter(proteinExpression, mVSamples, c = tissues)
        ax[index].set_xlabel(protein)
        ax[index].set_ylabel('Mean of Sample')
        ax[index].set_title(f"Protein expression across samples, by Tissue, $r∈{limitsR}$, $β_{{GLS}}∈{limitsGLS}${sortingSymbol}, corum={corum}, lung={lung}")
        ax[index].tick_params(labelsize=16)
        legend = [plt.Line2D([0], [0], marker='o', color=TISSUEPALETTE[key], label=key) for key in TISSUEPALETTE]
        ax[index].legend(handles = legend, fontsize=8, framealpha=0.2)

    plt.savefig('../images/scatterPxvsMVTruePositives.png')
"""
