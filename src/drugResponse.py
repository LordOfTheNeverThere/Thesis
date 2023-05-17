import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from resources import PATH, read, DrugResponseMatrix



if __name__ == "__main__":

    drugRes: DrugResponseMatrix = read(PATH + '/datasetsTese/drugResponse.pickle.gz')
    print(drugRes)
    # cellDeathThresh = range(3, 26)
    # allCapableDrugs = []

    # for n in cellDeathThresh:
    #     capableDrugs = drugRes.binrise(n,False).shape[0]
    #     allCapableDrugs.append(capableDrugs)
    
    # plt.bar(cellDeathThresh, allCapableDrugs)
    # plt.title("Drugs per Cell responsiveness threshold")
    # plt.xlabel("Cell responsiveness threshold")
    # plt.ylabel("Ammount of Drugs")

    # plt.savefig('../images/numOfCapableDrugsPerCellDeathThresh.png')´

    drugRes:pd.DataFrame = drugRes.binrise(0,False)
    sumCellDrugs = drugRes.sum(axis=1) # We get the sum of all cell lines that are killed by a drug (col)
    print(sumCellDrugs.sort_values(ascending=False))
    fig, ax = plt.subplots(1)
    plt.hist(sumCellDrugs.reset_index(drop=True), 20)
    ax.set_ylabel('# of Drugs')
    ax.set_xlabel('# Responsive Cell lines')
    ax.set_title("Histogram of Cell lines killed per drug")


    plt.savefig('../images/hitogramDrugsKillCelllines.png')