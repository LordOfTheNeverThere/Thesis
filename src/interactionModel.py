import pandas as pd
from resources import PATH, read, DrugResponseMatrix, ResiduesMatrix, ProteinsMatrix, PairwiseCorrMatrix
import statsmodels.api as sm
import statsmodels.formula.api as smf
from typing import Iterable
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.multitest import multipletests
import time as t

class DRInteractionPxModel():

    def __init__(self, ppis:Iterable[tuple[str,str]], proteomics:ProteinsMatrix, drugRes:DrugResponseMatrix, M:pd.DataFrame|pd.Series, fitIntercept=True, copy_X=True, standardisePx = True):

        self.ppis = ppis
        self.proteomics = proteomics.data
        self.drugRes = drugRes.data
        self.M = M
        self.fitIntercept = fitIntercept
        self.copy_X = copy_X
        self.standardisePx = standardisePx

    def getLinearModels(self, YName, XName, drugName) -> tuple[RegressionResultsWrapper, RegressionResultsWrapper, int]:
        """Get the Linear Models (Larger and Smaller) for the given Protein X and Y Names
        It does this by subsetting the proteomics, drugRes, and M dataframes to only include samples common to all dataframes.
        And then builds the OLS models from statsmodels.api library.

        Args:
            YName (str): Protein Y name
            XName (str): Protein X name
            drugName (str): Drug name

        Returns:
            sm.OLS: larger Linear Model from statsmodels.api library
            sm.OLS: Smaller Linear Model from statsmodels.api library
            int: number of samples common to all dataframes
        """
        
        Y = self.proteomics[YName]
        X = self.proteomics[XName]
        M = self.M
        drugRes = self.drugRes[drugName]
        #get samples common to all dataframes
        samplesCommon = list(set(Y.index) & set(X.index) & set(M.index) & set(drugRes.index)) # samples common to all dataframes
        #subset dataframes to common samples
        Y = Y.loc[samplesCommon]
        X = X.loc[samplesCommon]

        #number of samples in common, n
        n = len(samplesCommon)
        
        if self.standardisePx: # Zscore X if standardisePx is True
            X = (X - X.mean()) / X.std()

        M = M.loc[samplesCommon]
        drugRes = drugRes.loc[samplesCommon]

        commonDF = pd.concat([Y, X, M, drugRes], axis=1) #concatenate dataframes to build formula
        commonDF.columns = ['Py', 'Px'] + list(M.columns) + ['dR'] #rename columns to build formula

        #build formula containing all the M vars as independent variables, separately
        largeFormula = "Py ~ Px + " + " + ".join(list(M.columns)) + " + dR + dR:Px"
        smallForumla = "Py ~ Px + " + " + ".join(list(M.columns))

        largeModel = smf.ols(formula=largeFormula, data=commonDF).fit()
        smallModel = smf.ols(formula=smallForumla, data=commonDF).fit()
        
        print("This the the larger model's summary: \n")
        print(largeModel.summary())
        print("This the the smaller model's summary: \n")
        print(smallModel.summary())

        return largeModel, smallModel, n

    def fit(self)->pd.DataFrame:
        """Fit each Px and Py pair towards every drug in the drugRes dataframe.
            Calculate the Log likelihood p-value for the null hypothesis that the smaller model is correct, so the larger model does not add any covariate which is statistically significant.
            Or so to say the wilk's or likelihood ratio test.

        Returns:
            pd.DataFrame: The results of the fitting process, with the following columns: 
                Py, Px, drug, n, intercept, PxBeta, adherentBeta, semiAdherentBeta, suspensionBeta, unknownBeta, drugResBeta, interactionBeta, logLikePValue, llStatistic
        """        
        numOFPPIs = len(self.ppis)
        numOfDrugs = self.drugRes.shape[1]
        res = {
            'Py':[], 
            'Px':[], 
            'drug':[], 
            'n':[], 
            'intercept':[], 
            'PxBeta':[], 
            'adherentBeta':[], 
            'semiAdherentBeta':[], 
            'suspensionBeta':[], 
            'unknownBeta':[], 
            'drugResBeta':[], 
            'interactionBeta':[], 
            'logLikePValue':[],
            'llStatistic':[], 
            'fdr':[]}

        for ppi in self.ppis:
            for index,drugName in enumerate(list(self.drugRes.columns)):
                if index == 0:
                    start = t.time()
                

                YName = ppi[0]
                XName = ppi[1]
                largeModel, smallModel, n = self.getLinearModels(YName, XName, drugName)
                # P-value for the null hypothesis that the small model is correct, so the larger model does not add any covariate which is statistically significant
                llStatistic, logLikelihoodpValue,_ = largeModel.compare_lr_test(smallModel)

                res['Py'].append(YName)
                res['Px'].append(XName)
                res['drug'].append(drugName)
                res['n'].append(n)
                res['intercept'].append(largeModel.params['Intercept'])
                res['PxBeta'].append(largeModel.params['Px'])
                res['adherentBeta'].append(largeModel.params['Adherent'])
                res['semiAdherentBeta'].append(largeModel.params['SemiAdherent'])
                res['suspensionBeta'].append(largeModel.params['Suspension'])
                res['unknownBeta'].append(largeModel.params['Unknown'])
                res['drugResBeta'].append(largeModel.params['dR'])
                res['interactionBeta'].append(largeModel.params['dR:Px'])
                res['logLikePValue'].append(logLikelihoodpValue)
                res['llStatistic'].append(llStatistic)
                res['fdr'].append(multipletests(logLikelihoodpValue, method="fdr_bh")[1])

                if index == 0:
                    end = t.time()
                    print(f"The estimated time of completion if {numOFPPIs * numOfDrugs * (end-start) / 60} minutes")


        res = pd.DataFrame(res)
        return res
        



    #py ~ M + px + drugres + px:drugres

drugRes = read(PATH + '/internal/drugResponses/drugResponse.pickle.gz')
drugRes.data = drugRes.data.T
samplesheet = pd.read_csv(PATH + '/internal/samplesheet.csv', index_col=0)
ogProteomics: ProteinsMatrix = read(PATH + '/internal/proteomics/ogProteomics.pickle.gz')
# using only corum ppis that we were able to recall, with high confidence
# vaeGLSPairwise: PairwiseCorrMatrix = read(PATH + '/internal/pairwiseCorrs/VAE/glsPairCorr.pickle.gz')
# ppisOfInterest = set(vaeGLSPairwise.data.query("pValue < 0.001 & corum == 1").copy().index)
# ppisOfInterest = {(ppi.split(';')[0], ppi.split(';')[1]) for ppi in ppisOfInterest}

ppisOfInterest = [('MRPL38', 'MRPL45'), ('MRPL45', 'MRPL38')]
drugRes.data = drugRes.data.iloc[:,0:1]

M = pd.get_dummies(samplesheet['growth_properties'])
M = M.rename(columns={'Semi-Adherent': 'SemiAdherent'})

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


model = DRInteractionPxModel(ppisOfInterest, ogProteomics, drugRes, M, standardisePx=False)
model.fit()

