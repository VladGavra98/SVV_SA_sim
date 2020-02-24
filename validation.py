import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def deflectionAlongX(y,z,case):
    if case=="bending":
        Df = pd.read_csv('validation_processed/bending_processed.csv', header=0,index_col=0)
    elif case=="jam_bent":
        Df = pd.read_csv('validation_processed/jam_bent_processed.csv', header=0,index_col=0)
    elif case=="jam_straight":
        Df = pd.read_csv('validation_processed/jam_straight_processed.csv', header=0,index_col=0)
    else:
        print("Case not found")
        return

    yFiltered = Df.loc[Df['y'] == y]
    zFiltered = yFiltered.loc[yFiltered['z'] == z]
    xSort = zFiltered.sort_values(by='x')

    print(xSort['x'])
    plt.plot(xSort['x'],xSort['u2Loc1'])
    plt.show()


    print("zFiltered",zFiltered)
    print("zFiltered y",zFiltered['y'])
    print("zFiltered z",zFiltered['z'])

    return

deflectionAlongX(0,0,"bending")


def maxStressAlongX(case):
    if case == "bending":
        Df = pd.read_csv('validation_processed/bending_processed.csv', header=0,index_col=0)
    elif case == "jam_bent":
        Df = pd.read_csv('validation_processed/jam_bent_processed.csv', header=0,index_col=0)
    elif case == "jam_straight":
        Df = pd.read_csv('validation_processed/jam_straight_processed.csv', header=0,index_col=0)
    else:
        print("Case not found")
        return

    xArray = Df.x.unique()
    xArray.sort()
    maxStress = np.zeros(len(xArray))
    for idx,xval in enumerate(xArray):
        print("xval",xval)
        rows = Df.loc[Df['x'] == xval]
        print(rows)
        print(rows[['element','sMisesLoc1','sMisesLoc2']])
        avgMises = rows[['sMisesLoc1','sMisesLoc2']].mean(axis=1)
        print(avgMises)
        maxMises = avgMises.max()
        maxStress[idx] = maxMises

    plt.plot(xArray,maxStress)
    plt.show()
    print(xArray)
    return
maxStressAlongX("bending")

