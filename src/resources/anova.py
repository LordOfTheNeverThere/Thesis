import statsmodels.formula.api as smf
import pandas as pd
import numpy as np

class Anova():
    
    def __init__(self, data:pd.DataFrame, useStatsModels:bool):
        self.data = data
        self.useStatsModels = useStatsModels
        
        if useStatsModels:
            print('Warning:\n Using StatsModels for ANOVA, this is not recommended for large datasets, use useStatsModels=False for better performance')

    def fitOneWay(self, x:str, y:str):
        """Fits the analysis of variance model to the data.

        Args:
            x (str): Name of the independent categorical variable
            y (str): Name of the dependent variable, must be numerical
        """

        # assert that x is categorical, and y numerical
        data = self.data.copy()
        assert data[y].dtype == 'float64', 'y must be numerical'
        

        # fit the model
        #TODO: Make both options output the same value, the first outputs the difference in mean in relation to the first group, 
        # the second outputs the difference in mean in relation to the mean of all groups
        if self.useStatsModels:

            model = smf.ols(f'{y} ~ C({x})', data=data).fit()
            self.params = self.model.params

            return  self
        else:
            assert data[x].dtype == 'object', 'x must be categorical'
            #knowing that tau_i = mu_i - mu, where mu_i is the mean of the i-th group, and mu is the mean of all groups
            #we can calculate the mean of all groups, and then calculate the mean of each group, and then calculate the tau_i
            
            
            #calculate the mean of all groups
            globalMean = data[y].mean()
            groupMeans = data.groupby(x)[y].mean()
            #calculate the tau_i
            self.params = groupMeans - globalMean

            return self

