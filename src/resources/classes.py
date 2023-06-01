from __future__ import annotations #  postpone evaluation of annotations
import pandas as pd
from itertools import combinations
import numpy as np
import math
from sklearn.metrics import auc
import pickle
import gzip
from scipy.special import stdtr
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2
from resources import *
import scipy.special as special
from scipy import linalg
from sys import platform

def pearsonR(x, y, *, alternative='two-sided'):
    r"""
    Pearson correlation coefficient and p-value for testing non-correlation, 
    from Scipy, modified by Miguel Dinis de Sousa so as to increase pValue precion in Unix Systems.

    The Pearson correlation coefficient [1]_ measures the linear relationship
    between two datasets. Like other correlation
    coefficients, this one varies between -1 and +1 with 0 implying no
    correlation. Correlations of -1 or +1 imply an exact linear relationship.
    Positive correlations imply that as x increases, so does y. Negative
    correlations imply that as x increases, y decreases.

    This function also performs a test of the null hypothesis that the
    distributions underlying the samples are uncorrelated and normally
    distributed. (See Kowalski [3]_
    for a discussion of the effects of non-normality of the input on the
    distribution of the correlation coefficient.)
    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Pearson correlation at least as extreme
    as the one computed from these datasets.

    Parameters
    ----------
    x : (N,) array_like
        Input array.
    y : (N,) array_like
        Input array.
    alternative : {'two-sided', 'greater', 'less'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        The following options are available:

        * 'two-sided': the correlation is nonzero
        * 'less': the correlation is negative (less than zero)
        * 'greater':  the correlation is positive (greater than zero)

        .. versionadded:: 1.9.0

    Returns
    -------
    result : `~scipy.stats._result_classes.PearsonRResult`
        An object with the following attributes:

        statistic : float
            Pearson product-moment correlation coefficient.
        pvalue : float
            The p-value associated with the chosen alternative.

        The object has the following method:

        confidence_interval(confidence_level=0.95)
            This method computes the confidence interval of the correlation
            coefficient `statistic` for the given confidence level.
            The confidence interval is returned in a ``namedtuple`` with
            fields `low` and `high`.  See the Notes for more details.

    Warns
    -----
    `~scipy.stats.ConstantInputWarning`
        Raised if an input is a constant array.  The correlation coefficient
        is not defined in this case, so ``np.nan`` is returned.

    `~scipy.stats.NearConstantInputWarning`
        Raised if an input is "nearly" constant.  The array ``x`` is considered
        nearly constant if ``norm(x - mean(x)) < 1e-13 * abs(mean(x))``.
        Numerical errors in the calculation ``x - mean(x)`` in this case might
        result in an inaccurate calculation of r.

    See Also
    --------
    spearmanr : Spearman rank-order correlation coefficient.
    kendalltau : Kendall's tau, a correlation measure for ordinal data.

    Notes
    -----
    The correlation coefficient is calculated as follows:

    .. math::

        r = \frac{\sum (x - m_x) (y - m_y)}
                 {\sqrt{\sum (x - m_x)^2 \sum (y - m_y)^2}}

    where :math:`m_x` is the mean of the vector x and :math:`m_y` is
    the mean of the vector y.

    Under the assumption that x and y are drawn from
    independent normal distributions (so the population correlation coefficient
    is 0), the probability density function of the sample correlation
    coefficient r is ([1]_, [2]_):

    .. math::
        f(r) = \frac{{(1-r^2)}^{n/2-2}}{\mathrm{B}(\frac{1}{2},\frac{n}{2}-1)}

    where n is the number of samples, and B is the beta function.  This
    is sometimes referred to as the exact distribution of r.  This is
    the distribution that is used in `pearsonr` to compute the p-value.
    The distribution is a beta distribution on the interval [-1, 1],
    with equal shape parameters a = b = n/2 - 1.  In terms of SciPy's
    implementation of the beta distribution, the distribution of r is::

        dist = scipy.stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)

    The default p-value returned by `pearsonr` is a two-sided p-value. For a
    given sample with correlation coefficient r, the p-value is
    the probability that abs(r') of a random sample x' and y' drawn from
    the population with zero correlation would be greater than or equal
    to abs(r). In terms of the object ``dist`` shown above, the p-value
    for a given r and length n can be computed as::

        p = 2*dist.cdf(-abs(r))

    When n is 2, the above continuous distribution is not well-defined.
    One can interpret the limit of the beta distribution as the shape
    parameters a and b approach a = b = 0 as a discrete distribution with
    equal probability masses at r = 1 and r = -1.  More directly, one
    can observe that, given the data x = [x1, x2] and y = [y1, y2], and
    assuming x1 != x2 and y1 != y2, the only possible values for r are 1
    and -1.  Because abs(r') for any sample x' and y' with length 2 will
    be 1, the two-sided p-value for a sample of length 2 is always 1.

    For backwards compatibility, the object that is returned also behaves
    like a tuple of length two that holds the statistic and the p-value.

    References
    ----------
    .. [1] "Pearson correlation coefficient", Wikipedia,
           https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    .. [2] Student, "Probable error of a correlation coefficient",
           Biometrika, Volume 6, Issue 2-3, 1 September 1908, pp. 302-310.
    .. [3] C. J. Kowalski, "On the Effects of Non-Normality on the Distribution
           of the Sample Product-Moment Correlation Coefficient"
           Journal of the Royal Statistical Society. Series C (Applied
           Statistics), Vol. 21, No. 1 (1972), pp. 1-12.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> res = stats.pearsonr([1, 2, 3, 4, 5], [10, 9, 2.5, 6, 4])
    >>> res
    PearsonRResult(statistic=-0.7426106572325056, pvalue=0.15055580885344558)
    >>> res.confidence_interval()
    ConfidenceInterval(low=-0.9816918044786463, high=0.40501116769030976)

    There is a linear dependence between x and y if y = a + b*x + e, where
    a,b are constants and e is a random error term, assumed to be independent
    of x. For simplicity, assume that x is standard normal, a=0, b=1 and let
    e follow a normal distribution with mean zero and standard deviation s>0.

    >>> rng = np.random.default_rng()
    >>> s = 0.5
    >>> x = stats.norm.rvs(size=500, random_state=rng)
    >>> e = stats.norm.rvs(scale=s, size=500, random_state=rng)
    >>> y = x + e
    >>> stats.pearsonr(x, y).statistic
    0.9001942438244763

    This should be close to the exact value given by

    >>> 1/np.sqrt(1 + s**2)
    0.8944271909999159

    For s=0.5, we observe a high level of correlation. In general, a large
    variance of the noise reduces the correlation, while the correlation
    approaches one as the variance of the error goes to zero.

    It is important to keep in mind that no correlation does not imply
    independence unless (x, y) is jointly normal. Correlation can even be zero
    when there is a very simple dependence structure: if X follows a
    standard normal distribution, let y = abs(x). Note that the correlation
    between x and y is zero. Indeed, since the expectation of x is zero,
    cov(x, y) = E[x*y]. By definition, this equals E[x*abs(x)] which is zero
    by symmetry. The following lines of code illustrate this observation:

    >>> y = np.abs(x)
    >>> stats.pearsonr(x, y)
    PearsonRResult(statistic=-0.05444919272687482, pvalue=0.22422294836207743)

    A non-zero correlation coefficient can be misleading. For example, if X has
    a standard normal distribution, define y = x if x < 0 and y = 0 otherwise.
    A simple calculation shows that corr(x, y) = sqrt(2/Pi) = 0.797...,
    implying a high level of correlation:

    >>> y = np.where(x < 0, x, 0)
    >>> stats.pearsonr(x, y)
    PearsonRResult(statistic=0.861985781588, pvalue=4.813432002751103e-149)

    This is unintuitive since there is no dependence of x and y if x is larger
    than zero which happens in about half of the cases if we sample x and y.

    """
    n = len(x)
    if n != len(y):
        raise ValueError('x and y must have the same length.')

    if n < 2:
        raise ValueError('x and y must have length at least 2.')

    x = np.asarray(x)
    y = np.asarray(y)

    if (np.issubdtype(x.dtype, np.complexfloating)
            or np.issubdtype(y.dtype, np.complexfloating)):
        raise ValueError('This function does not support complex data')

    # If an input is constant, the correlation coefficient is not defined.
    if (x == x[0]).all() or (y == y[0]).all():
        msg = ("An input array is constant; the correlation coefficient "
               "is not defined.")
        print("An input array is constant; the correlation coefficient "
               "is not defined.")
        result = (np.nan, np.nan)
        return result

    # dtype is the data type for the calculations.  This expression ensures
    # that the data type is at least 64 bit floating point.  It might have
    # more precision if the input is, for example, np.longdouble.
    dtype = type(1.0 + x[0] + y[0])

    if n == 2:
        r = dtype(np.sign(x[1] - x[0])*np.sign(y[1] - y[0]))
        result = (r, 1.0)
        return result

    xmean = x.mean(dtype=dtype)
    ymean = y.mean(dtype=dtype)

    # By using `astype(dtype)`, we ensure that the intermediate calculations
    # use at least 64 bit floating point.
    xm = x.astype(dtype) - xmean
    ym = y.astype(dtype) - ymean

    # Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
    # scipy.linalg.norm(xm) does not overflow if xm is, for example,
    # [-5e210, 5e210, 3e200, -3e200]
    normxm = linalg.norm(xm)
    normym = linalg.norm(ym)

    threshold = 1e-13
    if normxm < threshold*abs(xmean) or normym < threshold*abs(ymean):
        # If all the values in x (likewise y) are very close to the mean,
        # the loss of precision that occurs in the subtraction xm = x - xmean
        # might result in large errors in r.
        print("An input array is nearly constant; the computed "
               "correlation coefficient may be inaccurate.")

    r = np.dot(xm/normxm, ym/normym)

    # Presumably, if abs(r) > 1, then it is only some small artifact of
    # floating point arithmetic.
    r = max(min(r, 1.0), -1.0)

    # As explained in the docstring, the p-value can be computed as
    #     p = 2*dist.cdf(-abs(r))
    # where dist is the beta distribution on [-1, 1] with shape parameters
    # a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
    # on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
    # shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
    # becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
    # to avoid a TypeError raised by btdtr when r is higher precision.)
    ab = n/2 - 1
    if platform == 'linux' or platform == 'darwin':

        if alternative == 'two-sided':
            prob = 2*special.btdtr(ab, ab, 0.5*(1 - abs(np.float128(r))))
        elif alternative == 'less':
            prob = 1 - special.btdtr(ab, ab, 0.5*(1 - abs(np.float128(r))))
        elif alternative == 'greater':
            prob = special.btdtr(ab, ab, 0.5*(1 - abs(np.float128(r))))
        else:
            raise ValueError('alternative must be one of '
                            '["two-sided", "less", "greater"]')
    else:
        
        if alternative == 'two-sided':
            prob = 2*special.btdtr(ab, ab, 0.5*(1 - abs(np.float64(r))))
        elif alternative == 'less':
            prob = 1 - special.btdtr(ab, ab, 0.5*(1 - abs(np.float64(r))))
        elif alternative == 'greater':
            prob = special.btdtr(ab, ab, 0.5*(1 - abs(np.float64(r))))
        else:
            raise ValueError('alternative must be one of '
                            '["two-sided", "less", "greater"]')

    return (r, prob)




class MatrixData:
    def __init__(self, filepath: str = None, data: pd.DataFrame = None, **readerKwargs):
        self.data = data
        assert filepath or (
            data is not None), 'There should be either a filepath or data'

        if filepath:
            self.data: pd.DataFrame = pd.read_csv(filepath, **readerKwargs)
            
        elif data is not None:
            self.data: pd.DataFrame = data.copy()
        
    def __str__(self) -> str:
        return str(self.data)

    def write(self, filepath:str):

        with gzip.open(filepath, 'wb') as f:
            pickle.dump(self,f)
        f.close()

    def query(self, query: str, inplace : bool = False) -> pd.DataFrame|None:


        if not inplace:
            return self.data.query(query).copy()
        else:
            self.data = self.data.query(query).copy()

    def compare(self, other: MatrixData, querySelf: str, queryOther: str, key: str = 'PPI') -> pd.DataFrame:
        """Query two PairwiseCorrMatrices with independent queries and get the merged dataframe as result

        Args:
            other (PairwiseCorrMatrix): another PairwiseCorrMatrix object
            querySelf (str): Query to apply on the self object
            queryOther (str): Query to apply on the other object
            key (str, optional): In what column name should the merge be happening on, acting as consesual index. Defaults to 'PPI'.

        Returns:
            pd.DataFrame: _description_
        """

        left: pd.DataFrame = self.query(querySelf).copy()
        right: pd.DataFrame = other.query(queryOther).copy()

        return left.merge(right, on=key)
class ppiDataset(MatrixData):

    def __init__(self, filepath:str = None, data: pd.DataFrame = None, proteinLabels: list = [], **readerKwargs):


        super().__init__(filepath, data, **readerKwargs)
        self.proteinLabels = proteinLabels
        self.ppis = set()


    def getPPIs(self, dataset:str) -> set:
        """Get the curated and observed ppis of a certains external PPI dataset

        Args:
            dataset (str): name of the dataset. It can either be 'corum', 'biogrid' or 'string'

        Returns:
            set: Set of tuples with proteinA and proteinB as putative PPI pairs
        """

        data = self.data.copy()
        ppiSet = None
        allowedDatasets = ['corum', 'biogrid', 'string']

        assert dataset in allowedDatasets, f"dataset not supported use one of the following 'corum', 'biogrid', 'string', got: {dataset}"

        if dataset == allowedDatasets[0]:

            def combinationsOfProteins(complx):

                if ';' in list(complx['subunits(Gene name)']):
                    complx['proteinTuple'] = list(combinations(
                        complx['subunits(Gene name)'].split(';'), 2))

                    return complx

            data = data.apply(
                lambda complx: combinationsOfProteins(complx), axis=1)
            ppiList = list(data.dropna()['proteinTuple'])
            ppiSet = {item for sublist in ppiList for item in sublist}

        elif (dataset == allowedDatasets[1]):

            # Filter Biogrid for certain parameters

            # Only allow ppis documented by physical interaction      
            # Filter out Homedymers which are not object of our study
            data = data.query("`Experimental System Type` == 'physical' and `Official Symbol Interactor A` != `Official Symbol Interactor B`").copy()

            data['proteinTuple'] = list(zip(data['Official Symbol Interactor A'], data['Official Symbol Interactor B']))
            ppiSet = set(data['proteinTuple'])

            

        elif (dataset == allowedDatasets[2]):
            data['proteinTuple'] = list(zip(data['proteinA'], data['proteinB']))
            ppiSet = set(data['proteinTuple'])

        
        self.ppis = ppiSet
        return ppiSet


class ProteinsMatrix(MatrixData):

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, **readerKwargs):

        super().__init__(filepath, data, **readerKwargs)

    def pearsonCorrelations(self, columnName: str, thresholdInteraction:int = 5) -> PairwiseCorrMatrix:
        """Calculate the pearson correlations and corresponding p-value, displaying them in a pairwise manner, returning an instance of the PairwiseCorrMatrix class

        Args:
            columnName (str): Name given to the df column with the correlation metric
            thresholdInteraction (int, optional): The minimum number of coincident samples for it to be considered a putative PPI

        Returns:
            PairwiseCorrMatrix: final data structure with all the information regarding pairwise correlations
        """
        if platform == 'linux' or platform == 'darwin':
            print('You are on an Unix System, using np.foat128 for greater precision on p-value')
        else:
            print('You need to be on an Unix System to take advantage of np.float128 for greater precision on p-value')

        data = self.data.copy()
        permutationsColumns = list(combinations(data.columns, 2))
        pairwiseCorr:dict = {'PPI':[],'proteinA': [],'proteinB': [], columnName: [], 'pValue':[], 'counts':[]}
        
        for (proteinAName, proteinBName) in permutationsColumns:

            proteinA = data[proteinAName].dropna(axis=0)
            proteinB = data[proteinBName].dropna(axis=0)
            samples  = proteinA.index.intersection(proteinB.index)
            proteinA = proteinA.loc[samples]
            proteinB = proteinB.loc[samples]

            count = len(proteinA)
            if count  < thresholdInteraction :
                continue
            
            
            (corr, pValue) = pearsonR(proteinA, proteinB)
            pairwiseCorr['PPI'].append(proteinAName + ';' + proteinBName)
            pairwiseCorr['proteinA'].append(proteinAName)
            pairwiseCorr['proteinB'].append(proteinBName)
            pairwiseCorr['pearsonR'].append(corr)
            pairwiseCorr['pValue'].append(pValue)
            pairwiseCorr['counts'].append(count)

        index = pairwiseCorr.pop('PPI')
        return PairwiseCorrMatrix(None, pd.DataFrame(pairwiseCorr, index=index))
      
    def tlsResidues(self, pairwiseCorr:PairwiseCorrMatrix) -> ResiduesMatrix:

        pairwiseCorr = pairwiseCorr.data
        proteomics = self.data.copy()
        tlsResList = []
        ppis =  pairwiseCorr.loc[pairwiseCorr['corum'] == 1].loc[pairwiseCorr['glsCoefficient'] > 0.8].index 
        #Get the ppis that are most likely true ppis, so that we can analyse what samples do not correspond to the correlation, 
        # are farthest from the linear regession line, hence, have greatest TLS and so are samples of interest where the PPI likely. 
        # So it would be interesting to se afterwards if that sample has a responsiveness to a drug all the other samples do not meaning 
        # we are in a presence of a PPI that might be correlated to a feature, a certain drug responsiveness
        

        for ppi in ppis:
            proteinA = ppi.split(';')[0]
            proteinB = ppi.split(';')[1]

            X = proteomics.loc[:,proteinA].dropna(axis=0) #Get X and Y data
            Y = proteomics.loc[:,proteinB].dropna(axis=0)
            samplesInCommon = X.index.intersection(Y.index)

            if len(samplesInCommon) < 5: # We have no important information from protein protein interactions with less than 5 coocorence
                
                print(f"{ppi} was not used for calculating the tls Redidues Matrix because it did not have at least 5 samples")
                continue

            X=X.loc[samplesInCommon] #Locate samples that both have some value (not nan)
            Y=Y.loc[samplesInCommon]

            meanX = X.mean()
            meanY = Y.mean()

            meanErrorX = X - meanX #Start Calculating quantities to minimise for the tlsCoef Calculation
            meanErrorY = Y - meanY

            meanSqErrorX = meanErrorX ** 2
            meanSqErrorY = meanErrorY ** 2

            u = meanSqErrorX.sum()
            v = meanSqErrorY.sum()
            r = (meanErrorX * meanErrorY).sum()
            w = v - u

            tlsCoef = (w + (w**2 + r**2)**0.5) #Calculating tls Coefficient
            tlsCoef = tlsCoef/r

            intercept = meanY - (tlsCoef * meanX) #Intercept of linear fit
            predY = intercept + (tlsCoef * X)
            residues = Y - predY # TLS Residues in absolute val
            residues = pd.DataFrame(residues,columns=[ppi])
            tlsResList.append(residues)


        tlsResData = pd.concat(tlsResList, join='outer', sort=False)
        return ResiduesMatrix(None,tlsResData)

    def getGLSCorr(self, pValues: bool = True, listCovMatrix:list[pd.DataFrame] = None, coefColumnName :str = 'glsCoefficient') -> PairwiseCorrMatrix:
        """Get the GLS coeficents between each Protein X and Y, where X != Y, these will measure the correlation between each protein. 
        But this method differs from the pearsonCorrelations since it has underneath a GLM where the covariance matrix can be any specified.
        This covariance matrix will transform both X and y of the proteinData.data. By default this covariance matrix is calculated with proteinData.data
        Where we get the covariance between samples, as a similarity measure between samples. This tweak is speacilly important if our residuals
        correlate with X, meaning we are not in the normal OLS case.

        Args:
            pValues (bool, optional): Add the pValues of each gls Coefficient to the output data. Defaults to True.
            listCovMatrix (list[pd.DataFrame], optional):List of matrices to use to calculate covariance, if only one is required insert [matrix]. Defaults to None.
            coefColumnName (str, optional): Name to appear on the Dataframe' Column of glsCoefs. Defaults to 'glsCoefficient'.

        Returns:
            PairwiseCorrMatrix: Data structure with all above information
        """
        proteinData = self.data.copy()

        if listCovMatrix is not None:
            # The matrix used for the covariance is not the one used as an X in the linear regression, 
            # we used this to allow the use of the genomic matrix as the matrix where the cov of each sample would be calculated, 
            # since the tecidual bias of proteomic expression would be more present in the genome, 
            # since this layer is upstream and more connected to the information the cell has from its diferentiation, 
            # so it becomes easier to distinguish samples by tissue, by using its genomic data. 
            # Therefore calculating the cov of genomics.csv could give a more correct value on the true covariation of two samples 
            # and correct the possible correlation of PPIs that are simply due to the diferent baseline expression of proteins belonging in the ppi

            listCovMatrix.append(proteinData)
            samplesInCommon = pd.concat(listCovMatrix, axis=1, join='inner').index #Conactenate all matrices and select their indices which would be the common ones
            listCovMatrix.pop(-1) # Remove recently added proteinData to the list

            proteinData = proteinData.loc[samplesInCommon, :] 
            proteinDataMean = proteinData.fillna(proteinData.mean())
            proteinDataMean.dropna(axis=1, inplace=True) #We delete all columns with nan values because after all this preprocessing they must be completetly empty columns
            covMatrix = np.zeros((len(samplesInCommon), len(samplesInCommon)))

            for dataForCov in listCovMatrix:
                
                dataForCov = dataForCov.loc[samplesInCommon, :]
                dataForCov.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #How we are handling missing data, there should be at leats 20% of missingness for a collumn to be dropable
                dataForCov = dataForCov.fillna(dataForCov.mean())
                # calculate covariance matrix in order to see the covariace between samples, and notice tecidual patterns codified in the samples
                dataForCov = np.cov(dataForCov)
                covMatrix = covMatrix + dataForCov




        else: # The matrix used for the covariance is the same as that used for X in linear regression, we are whitening while taking into account the covariation of our data in X
            
            proteinData.dropna(axis=1, thresh=round(proteinData.shape[0] * 0.2), inplace=True) #We require that a protein has about 20% missingness for it to be considered a dropable column
            proteinDataMean = proteinData.fillna(proteinData.mean())
            dataForCov = proteinDataMean
            # calculate covariance matrix in order to see the covariace between samples, and notice tecidual patterns codified in the samples
            covMatrix = np.cov(dataForCov)


        proteinNames = proteinDataMean.columns.str.split(' ').str.get(0).to_numpy()
        proteinNames = [protein1 + ';' + protein2 for i, protein1 in enumerate(proteinNames)  for j, protein2 in enumerate(proteinNames) if j > i]


        #invert it the covariance matrix
        covMatrix = np.linalg.inv(covMatrix)

        # Decompose it with Cholesky, returning the lower triangular matrix of the positive definite matrix covMatrix, because cov(x1,x2) == cov(x2,x1)
        cholsigmainvMean = np.linalg.cholesky(covMatrix)


        # Whittening transformation, we codify our data into a space where each the variance of each covariate is the same and equal to one, 
        # so we are kind like normalising it, in fact that's exactly what we are doing ~ N(0,I) As they call it warping...
        warpedProteinsMean = proteinDataMean.T.values @ cholsigmainvMean

        # The Mean of the variances of a sample across all gene symbols is used as a beta0 for linear regression
        warpedIntereceptMean = cholsigmainvMean.T.sum(axis=0)


        def linear_regression(warped_screens, warped_intercept):
            GLS_coef = np.empty((len(warped_screens), len(warped_screens)))
            GLS_se = np.empty((len(warped_screens), len(warped_screens)))
            ys = warped_screens.T
            for proteinIndex in range(len(warped_screens)):
                
                X = np.stack((warped_intercept, warped_screens[proteinIndex]), axis=1)
                if np.any(np.isnan(X)) :
                    print(proteinIndex)
                    print('\n')
                    print(X)
                coef, residues = np.linalg.lstsq(X, ys, rcond=None)[:2]
                df = warped_screens.shape[1] - 2
                GLS_coef[proteinIndex] = coef[1]
                GLS_se[proteinIndex] = \
                    np.sqrt(np.linalg.pinv(X.T @ X)[1, 1] * residues / df)
            return GLS_coef, GLS_se

        GLS_coef, GLS_se = linear_regression(warpedProteinsMean, warpedIntereceptMean)

        df = warpedProteinsMean.shape[1] - 2

        #   Construct new PairwiseGLSCoefs Matrix


        glsCoefs = GLS_coef[np.triu_indices(GLS_coef.shape[0], k=1)]
        if pValues: #We might not want to add pValues so that we use less memory
            GLS_p = 2 * stdtr(df, -np.abs(GLS_coef / GLS_se))
            np.fill_diagonal(GLS_p, 1)
            glsPValues = GLS_p[np.triu_indices(GLS_p.shape[0], k=1)]
            pairwiseCorrData = pd.DataFrame(
                {coefColumnName: glsCoefs, 'pValue': glsPValues}, index=proteinNames)
        else:
            pairwiseCorrData = pd.DataFrame(
                {coefColumnName: glsCoefs}, index=proteinNames)
        
        pairwiseCorrData = PairwiseCorrMatrix(None, pairwiseCorrData)
        pairwiseCorrData.data.index.name = 'PPI'

        return pairwiseCorrData


class PairwiseCorrMatrix(MatrixData):

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, aucs = dict(), ** readerKwargs):
        """_summary_

        Args:
            filepath (str, optional): filepath where Datafrane is stored to instatiate the object. Defaults to None.
            data (pd.DataFrame, optional):Dataframe to instatiate Object. Defaults to None.
            aucs (_type_, optional): dict with name of columns to which we calculate the auc and value. Defaults to None.
        """
        super().__init__(filepath, data, **readerKwargs)

        self.corrCumSums = aucs
        self.indexes = aucs
        self.aucs:dict = aucs 
        self.labels = aucs

    def __str__(self):

        print = super().__str__()
        for columnName, aucVal in self.aucs.items():

            print = print + '\n' +str(aucVal) + ' ' +str(columnName) + '\n' + str(self.labels[columnName])
        

        return print

    def addGroundTruth(self, ppis: set, externalDatasetName: str):
        """Append the binary values of a putative PPI, from an external dataset (e.g Corum), to our pairwise correlation Dataframe

        Args:
            ppis (ppiDataset): ppiDataset of the external ppi dataset used
            data (pd.DataFrame): Pairwise correlation dataframe
            externalDatasetName (str): Name to give to the binary column holding the truth value of an PPI is seen in that external Dataset
            filepath (str): The name to give to the final file, with the added binary column

        Returns:
            _type_: Data with added column
        """
        data = self.data.copy()
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


        self.data = data

        return data
    
    def aucCalculator(self, yColumnName:str, label:str, proxyColumn:str, ascending:bool ):
        """Adds the value of AUC of the Recall curve using a specified external PPI dataset with yColumnName

        Args:
            yColumnName (str): Name of the df column where there is the truth value of the existence or not of the PPI in the reported external PPI dataset
            label (str): Text which will show up as label next to the value of the AUC, e.g 'Baseline Auc == 0.9' 
            proxyColumn (str): Name of the column of the statistical meausure to quantify the probability of PPI existence
            ascending(bool): should the proxyColumn be ordered ascending or not for best AUC calculation
        """
        pairwiseCorr = self.data 

        
        pairwiseCorr.sort_values(by=proxyColumn, ascending=ascending, inplace=True) # We sort rows by the smallest to greatest pValues

        self.corrCumSums[proxyColumn] = np.cumsum(
            pairwiseCorr[yColumnName]) / np.sum(pairwiseCorr[yColumnName])
        
        self.indexes[proxyColumn] = np.array(pairwiseCorr.reset_index().index) / \
            pairwiseCorr.shape[0]
        self.aucs[proxyColumn] = auc(self.indexes[proxyColumn], self.corrCumSums[proxyColumn]) # update aucs dict to have a new auc for a specific proxy column

        if not label: #if the user did not insert any label default it
            self.labels[proxyColumn] = f"(AUC {proxyColumn} {self.aucs[proxyColumn]:.2f})"
        
        self.labels[proxyColumn] =  label + f" (AUC {proxyColumn} {self.aucs[proxyColumn]:.2f})"

    def aucsCalculator(self, yColumnNameList:list[str], label:str, proxyColumnList:list[str], ascendingList:list[bool], filepath:str = None ):

        for aucIndex in range(len(yColumnNameList)):
            self.aucCalculator(yColumnNameList[aucIndex],label,proxyColumnList[aucIndex], ascendingList[aucIndex])
        
    
        if filepath is not None:
            self.write(filepath)
        print(self)

    
class DrugResponseMatrix(MatrixData):
    """Class interface and methods for the drug response data"""

    def __init__(self, filepath: str=None, data: pd.DataFrame=None, **readerKwargs):
        super().__init__(filepath, data, **readerKwargs)


    def binrise(self, deathThresh: int = 3, inplace: bool = True):
        """Creates a Binary representation of the Drug response matrix, where 0 means no efficacy and 1 efficacy. This is done bye calculating a threshold, 50% of the natural log of the [max screen] 


        Args:
            deathThresh (int, optional): Number required of 'responsive' cell lines for a certain drug to be considered in the final set. Defaults to 3.
        """
        data = self.data.copy()
        maxScreenConc = pd.read_csv(PATH + '/externalDatasets/drugMaxScreenConcentration.csv', index_col='Unnamed: 0')
        maxScreenConc.index.name = 'drug'
        
        data = data.merge(maxScreenConc, on='drug') # We only work with drugs for which we have a max screen concentration and IC50
        # With 50% of the natural log of the max Screen Concentration we can create an empirically proven threshold, 
        # where generelly IC50 values greater than that mean that The drug has no efficacy, and below drug has efficancy, 
        # or in other words the cell line doesn't or does 'respond to a drug'
        data['efficacyThreshold']= data['MAX_CONC'].apply(lambda row: math.log(row * 0.5)) 
        data.drop(columns=['MAX_CONC'], inplace=True)


        # Create mask that see if all columns are less than the values in the threshold col, and convert them to ints so the bool's become either 0 or 1
        data:pd.DataFrame = data.apply(lambda row: row < row['efficacyThreshold'], axis=1).astype(int)

        relevantDrugs = (data.sum(axis=1) >= deathThresh) # Condition that only accounts for drugs that kill at leats 3 cell lines

        if inplace:
            self.data = data.loc[relevantDrugs].drop(columns=['efficacyThreshold']) # apply condition (We lose 130 drugs, 697 -> 567) 26/4/23
        else:
            return data.loc[relevantDrugs].drop(columns=['efficacyThreshold'])
        
class ResiduesMatrix(MatrixData):

    def __init__(self, filepath: str=None, data: pd.DataFrame=None, **readerKwargs):
        super().__init__(filepath, data, **readerKwargs)



class LinearRegression:
    """Authored by professor Emanuel Gonçalves, this call allows fo the computation of the linear regression coefs for a set of features. 
    Do check if the response variable is Normally distributed, since this is nor a GLM
    """
    
    

    def __init__(
        self,
        Y,
        X,
        M,
        M2=None,
        normalize=False,
        fit_intercept=True,
        copy_X=True,
        n_jobs=4,
        verbose=1,
    ):
        self.samples = set.intersection(
            set(Y.index),
            set(X.index),
            set(M.index),
            set(Y.index) if M2 is None else set(M2.index),
        )

        self.X = X.loc[self.samples]
        self.X = self.X.loc[:, self.X.count() > (M.shape[1] + (1 if M2 is None else 2))]
        self.X_ma = np.ma.masked_invalid(self.X.values)

        self.Y = Y.loc[self.samples]
        self.Y = self.Y.loc[:, self.Y.std() > 0]

        self.M = M.loc[self.samples]

        self.M2 = M2.loc[self.samples, self.X.columns] if M2 is not None else M2

        self.normalize = normalize
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.n_jobs = n_jobs

        self.verbose = verbose
        # self.log = logging.getLogger("Crispy")

    def model_regressor(self):
        regressor = LinearRegression(
            fit_intercept=self.fit_intercept,
            normalize=self.normalize,
            copy_X=self.copy_X,
            n_jobs=self.n_jobs,
        )
        return regressor

    @staticmethod
    def loglike(y_true, y_pred):
        nobs = len(y_true)
        nobs2 = nobs / 2.0

        ssr = np.power(y_true - y_pred, 2).sum()

        llf = -nobs2 * np.log(2 * np.pi) - nobs2 * np.log(ssr / nobs) - nobs2

        return llf

    @staticmethod
    def multipletests_per(
        associations, method="fdr_bh", field="pval", fdr_field="fdr", index_cols=None
    ):
        index_cols = ["y_id"] if index_cols is None else index_cols

        d_unique = {tuple(i) for i in associations[index_cols].values}

        df = associations.set_index(index_cols)

        df = pd.concat(
            [
                df.loc[i]
                .assign(fdr=multipletests(df.loc[i, field], method=method)[1])
                .rename(columns={"fdr": fdr_field})
                for i in d_unique
            ]
        ).reset_index()

        return df

    def fit_matrix(self):
        lms = []

        for x_idx, x_var in enumerate(self.X):
            if self.verbose > 0:
                self.log.info(f"LM={x_var} ({x_idx})")

            # Mask NaNs
            x_ma = np.ma.mask_rowcols(self.X_ma[:, [x_idx]], axis=0)

            # Build matrices
            x = self.X.iloc[~x_ma.mask.any(axis=1), [x_idx]]
            y = self.Y.iloc[~x_ma.mask.any(axis=1), :]

            # Covariate matrix (remove invariable features and add noise)
            m = self.M.iloc[~x_ma.mask.any(axis=1), :]
            if self.M2 is not None:
                m2 = self.M2.iloc[~x_ma.mask.any(axis=1), [x_idx]]
                m = pd.concat([m2, m], axis=1)
            m = m.loc[:, m.std() > 0]
            m += np.random.normal(0, 1e-4, m.shape)

            # Fit covariate model
            lm_small = self.model_regressor().fit(m, y)
            lm_small_ll = self.loglike(y, lm_small.predict(m))

            # Fit full model: covariates + feature
            lm_full_x = np.concatenate([m, x], axis=1)
            lm_full = self.model_regressor().fit(lm_full_x, y)
            lm_full_ll = self.loglike(y, lm_full.predict(lm_full_x))

            # Log-ratio test
            lr = 2 * (lm_full_ll - lm_small_ll)
            lr_pval = chi2(1).sf(lr)

            # Assemble + append results
            res = pd.DataFrame(
                dict(
                    y_id=y.columns,
                    x_id=x_var,
                    n=y.attrs["nan_mask"].loc[y.columns, x.index].sum(1) if "nan_mask" in y.attrs else len(x),
                    beta=lm_full.coef_[:, -1],
                    lr=lr.values,
                    covs=m.shape[1],
                    pval=lr_pval,
                    fdr=multipletests(lr_pval, method="fdr_bh")[1],
                )
            )

            lms.append(res)

        lms = pd.concat(lms, ignore_index=True).sort_values("pval")

        return lms

    @staticmethod
    def lm_residuals(y, x, fit_intercept=True, add_intercept=False):
        # Prepare input matrices
        ys = y.dropna()

        xs = x.loc[ys.index].dropna()
        xs = xs.loc[:, xs.std() > 0]

        ys = ys.loc[xs.index]

        if ys.shape[0] <= xs.shape[1]:
            return None

        # Linear regression models
        lm = LinearRegression(fit_intercept=fit_intercept).fit(xs, ys)

        # Calculate residuals
        residuals = ys - lm.predict(xs) - lm.intercept_

        # Add intercept
        if add_intercept:
            residuals += lm.intercept_

        return residuals
    






        

