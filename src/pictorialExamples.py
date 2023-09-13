import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import stdtr
import matplotlib.pyplot as plt
import matplotlib
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from resources import ProteinsMatrix, PATH, drawRecallCurves, read, ppiDataset, PairwiseCorrMatrix


if __name__ == '__main__':


    def plotPictorial(data:pd.DataFrame, filepath:str|None = None, xName:str = 'Px', yName:str = 'Py', DRName:str = 'DR'):

        # sort drug response so that the smallest shows up first, ascending==True
        data = data.sort_values(DRName, ascending=False)

        plottingData = data.loc[:,[xName, yName]]
        #standardize the data 
        plottingData = pd.DataFrame(StandardScaler().fit_transform(plottingData) ,columns=plottingData.columns, index=plottingData.index)

        #Add interactor
        plottingData = plottingData.join(data.loc[:,DRName], how = 'inner')


        plt.figure(figsize=(10, 10))
        scatter = sns.scatterplot(data=plottingData, x=xName, y=yName, hue=DRName, size=DRName, palette="viridis", alpha=1, edgecolor='none', s=15)
        norm = matplotlib.colors.Normalize(vmin=data[DRName].min(), vmax=data[DRName].max())
        # Add Colour Map
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        scatter.get_legend().remove()
        scatter.figure.colorbar(sm, label='Drug Response')

        plt.title(f'Protein expression \n with Drug Response')
        plt.xlabel(str(xName))
        plt.ylabel(str(yName))
        plt.tight_layout()

        if filepath is not None:
            plt.savefig(filepath)
        else:
            plt.show()
        plt.close()





    def getExample(
            proteinsCorrelated:bool, 
            normallyHighDR:bool, 
            clusterAbove:bool, 
            lineSamplesNum:int = 100, 
            clusterSamplesNum:int = 20):
        np.random.seed(0)
        
        results ={'index':[], 'Px':[], 'Py':[], 'DR':[]}

        if proteinsCorrelated:
            Px = np.random.normal(0, 1, lineSamplesNum)
            Px = Px + np.random.normal(0, 0.1, lineSamplesNum)
            Py = 1.5*Px + np.random.normal(0, 0.7, lineSamplesNum)

            PxCluster = np.random.normal(4, 1, clusterSamplesNum)
            PxCluster = PxCluster + np.random.normal(0, 0.1, clusterSamplesNum)

            if clusterAbove:
                
                PyCluster = 1.5*PxCluster + np.random.normal(0, 1, clusterSamplesNum) + 4

            else:

                PyCluster = 1.5*PxCluster + np.random.normal(0, 1, clusterSamplesNum) - 8
        else:
            Px = np.random.normal(0, 1, lineSamplesNum)
            Px = Px + np.random.normal(0, 0.1, lineSamplesNum)
            Py = -1.5*Px - np.random.normal(0, 0.7, lineSamplesNum)

            PxCluster = np.random.normal(4, 1, clusterSamplesNum)
            PxCluster = PxCluster + np.random.normal(0, 0.1, clusterSamplesNum)

            if clusterAbove:
                
                PyCluster = -1.5*PxCluster + np.random.normal(0, 1, clusterSamplesNum) + 8

            else:

                PyCluster = -1.5*PxCluster + np.random.normal(0, 1, clusterSamplesNum) - 4

        Px = np.concatenate((Px, PxCluster))
        Py = np.concatenate((Py, PyCluster))
        
        if normallyHighDR:
            DR = np.random.normal(10, 1, lineSamplesNum); DRCluster =  np.random.normal(5, 1, clusterSamplesNum)
        else:
            DR = np.random.normal(5, 1, lineSamplesNum); DRCluster =  np.random.normal(10, 1, clusterSamplesNum)

        DR = np.concatenate((DR, DRCluster))

        results['index'] = list(range(lineSamplesNum + clusterSamplesNum))
        results['Px'] = list(Px)
        results['Py'] = list(Py)
        results['DR'] = list(DR)


        return results
    

    for proteinsCorrelated in [True, False]:
        for normallyHighDR in [True, False]:
            for clusterAbove in [True, False]:
                results = getExample(proteinsCorrelated, normallyHighDR, clusterAbove)
                index = results.pop('index')
                results = pd.DataFrame(results, index=index)
                plotPictorial(results, f'pictorialExampleProteinsCorrelated={proteinsCorrelated}NormallyHigh={normallyHighDR}ClusterAbove={clusterAbove}.png')
