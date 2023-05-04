import pandas as pd
import utils
from classes import DrugResponseMatrix
import matplotlib.pyplot as plt
from env import PATH




if __name__ == "__main__":

    drugRes : DrugResponseMatrix = utils.read(PATH + '/datasetsTese/drugResponse.pickle.gz')
    cellDeathThresh = range(3, 26)
    allCapableDrugs = []

    for n in cellDeathThresh:
        capableDrugs = drugRes.binrise(n,False).shape[0]
        allCapableDrugs.append(capableDrugs)
    
    plt.bar(cellDeathThresh, allCapableDrugs)
    plt.title("Drugs per Cell responsiveness threshold")
    plt.xlabel("Cell responsiveness threshold")
    plt.ylabel("Ammount of Drugs")

    plt.savefig('../images/numOfCapableDrugsPerCellDeathThresh.png')