import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
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

    fig, ax = plt.subplots()
    plt.title("Aileron deflection VS distance in x-direction")
    plt.plot(xSort['x'], xSort['u1Loc1'], label="deflection in x")
    plt.plot(xSort['x'], xSort['u2Loc1'],label="deflection in y")
    plt.plot(xSort['x'], xSort['u3Loc1'], label="deflection in z")
    plt.xlabel("x [mm]")
    plt.ylabel("Deflection [mm]")
    plt.grid()
    plt.legend()
    ax.autoscale()
    plt.show()
    return

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
        rows = Df.loc[Df['x'] == xval]
        avgMises = rows[['sMisesLoc1','sMisesLoc2']].mean(axis=1)
        maxMises = avgMises.max()
        maxStress[idx] = maxMises

    fig, ax = plt.subplots()
    plt.title("Maximum Von Mises stress VS distance in x-direction")
    plt.plot(xArray,maxStress,label="maximum stress")
    plt.xlabel("x [mm]")
    plt.ylabel("Von Mises stress [Pa]")
    plt.grid()
    plt.legend()
    ax.autoscale()
    plt.show()
    return

def stressCrossSection(x,case):
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
    absolute_difference_function = lambda list_value: abs(list_value - x)
    nearestX = min(xArray, key=absolute_difference_function)
    print("Nearest x found is:",nearestX)
    rows = Df.loc[Df['x'] == nearestX]
    rows = rows.assign(avgMises=rows.loc[:,['sMisesLoc1', 'sMisesLoc2']].mean(axis=1))
    minAvgMises = min(rows.loc[:,'avgMises'])
    maxAvgMises = max(rows.loc[:,'avgMises'])
    diff = maxAvgMises-minAvgMises
    rows = rows.assign(c=((rows.loc[:,'avgMises']-minAvgMises)/diff))

    ### section 1
    rows1 = rows.loc[(rows['y']>=0) & (rows['z']<0)]
    rows1 = rows1.sort_values(by='y')
    print(rows1)
    #rows1.to_csv(r'test.csv')
    avg1 = rows1.groupby(['y','z']).mean().reset_index()
    #avg.to_csv(r'testavg.csv')
    #print("avg",avg)
    print(avg1['c'].to_numpy())


    multiline([avg1['z'].to_numpy()],[avg1['y'].to_numpy()],avg1['c'].to_numpy(),cmap='jet',lw=2)
    #plt.scatter(rows['z'],rows['y'],c=rows['avgMises'],cmap="jet")
    #minAvgMises = min(rows['avgMises'])



    return



def multiline(xs, ys, c, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    fig, ax = plt.subplots()
    plt.title("Cross section")
    # create LineCollection
    points = np.array([xs, ys]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='jet')
    # Set the values used for colormapping
    lc.set_array(c)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    fig.colorbar(line)


    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    plt.grid()
    plt.axis('equal')
    ax.autoscale()
    plt.show()
    return


#deflectionAlongX(0,0,"bending")
#maxStressAlongX("bending")
stressCrossSection(1074,"bending")
