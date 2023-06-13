import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from resources import PATH, read, DrugResponseMatrix



if __name__ == "__main__":

    drugRes: DrugResponseMatrix = read(PATH + '/internal/drugResponse.pickle.gz')
    # drugRes.data.plot.hist()
    # plt.savefig('../images/histOfDrugResponseDF.png')
