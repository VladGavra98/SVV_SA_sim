import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from scipy import spatial

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

    uniqueYZ = Df.drop_duplicates(subset=['y','z'])
    A = np.zeros((uniqueYZ['z'].size,2))
    A[:,0] = uniqueYZ['y']
    A[:, 1] = uniqueYZ['z']
    nearestYZ = A[spatial.KDTree(A).query([y,z])[1]]
    print("Nearest (y,z) found is:", nearestYZ)
    nearestY = nearestYZ[0]
    nearestZ = nearestYZ[1]
    print("Plotting deflection in (x,y,z)-direction for case: %s" %(case))


    yFiltered = Df.loc[Df['y'] == nearestY]
    zFiltered = yFiltered.loc[yFiltered['z'] == nearestZ]
    xSort = zFiltered.sort_values(by='x')

    fig, ax = plt.subplots()
    plt.suptitle("Aileron deflection VS distance in x-direction for case: %s"%(case),fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    ax.set_ylabel('y [mm]', fontsize=12.0)
    ax.set_xlabel('z [mm]', fontsize=12.0)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)

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

    print("Plotting maximum Von Mises stress in x-direction for case: %s" %(case))

    xArray = Df.x.unique()
    xArray.sort()
    maxStress = np.zeros(len(xArray))
    for idx,xval in enumerate(xArray):
        rows = Df.loc[Df['x'] == xval]
        avgMises = rows[['sMisesLoc1','sMisesLoc2']].mean(axis=1)
        maxMises = avgMises.max()
        maxStress[idx] = maxMises

    fig, ax = plt.subplots()
    plt.suptitle("Maximum Von Mises stress VS distance in x-direction for case: %s" % (case),fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    ax.set_ylabel('y [mm]', fontsize=12.0)
    ax.set_xlabel('z [mm]', fontsize=12.0)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)

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

    fig,ax = plt.subplots()

    xArray = Df.x.unique()
    #print(len(xArray))
    absolute_difference_function = lambda list_value: abs(list_value - x)
    nearestX = min(xArray, key=absolute_difference_function)
    print("Nearest x found is:",nearestX)
    print("Plotting Von Mises stress in cross section for case: %s" % (case))
    rows = Df.loc[Df['x'] == nearestX]
    rows = rows.assign(avgMises=rows.loc[:,['sMisesLoc1', 'sMisesLoc2']].mean(axis=1))
    minAvgMises = min(rows.loc[:,'avgMises'])
    maxAvgMises = max(rows.loc[:,'avgMises'])
    #diff = maxAvgMises-minAvgMises
    #rows = rows.assign(c=((rows.loc[:,'avgMises']-minAvgMises)/diff))
    norm = plt.Normalize(rows['avgMises'].min(), rows['avgMises'].max())

    #plotting paramters
    plt.suptitle("Von Mises stress [Pa] at cross section x=%d [mm] for case: %s" % (nearestX,case), fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    ax.set_ylabel('y [mm]', fontsize=12.0)
    ax.set_xlabel('z [mm]', fontsize=12.0)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    ### section 1
    rows1 = rows.loc[(rows['y']>=0) & (rows['z']<0)]
    rows1 = rows1.append(rows.loc[(rows['y']==102.5) & (rows['z']==0)])
    rows1 = rows1.sort_values(by='y')
    avg1 = rows1.groupby(['y','z']).mean().reset_index()
    points1 = np.array([avg1['z'].to_numpy(),avg1['y'].to_numpy()]).T.reshape(-1, 1, 2)
    segments1 = np.concatenate([points1[:-1], points1[1:]], axis=1)
    lc1 = LineCollection(segments1, cmap='jet', norm=norm)
    lc1.set_array(avg1['avgMises'])
    lc1.set_linewidth(3)
    ax.add_collection(lc1)

    ### section 2
    rows2 = rows.loc[(rows['z'] == 0)]
    rows2 = rows2.sort_values(by='y',ascending=False)
    avg2 = rows2.groupby(['y', 'z']).mean().reset_index()
    points2 = np.array([avg2['z'].to_numpy(), avg2['y'].to_numpy()]).T.reshape(-1, 1, 2)
    segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)
    lc2 = LineCollection(segments2, cmap='jet', norm=norm)
    lc2.set_array(avg2['avgMises'])
    lc2.set_linewidth(3)
    ax.add_collection(lc2)

    ### section 3
    rows3 = rows.loc[(rows['y'] <= 0) & (rows['z'] < 0)]
    rows3 = rows3.append(rows.loc[(rows['y'] == -102.5) & (rows['z'] == 0)])
    rows3 = rows3.sort_values(by='z')
    avg3 = rows3.groupby(['y', 'z']).mean().reset_index()
    points3 = np.array([avg3['z'].to_numpy(), avg3['y'].to_numpy()]).T.reshape(-1, 1, 2)
    segments3 = np.concatenate([points3[:-1], points3[1:]], axis=1)
    lc3 = LineCollection(segments3, cmap='jet', norm=norm)
    lc3.set_array(avg3['avgMises'])
    lc3.set_linewidth(3)
    ax.add_collection(lc3)

    ### section 4
    rows4 = rows.loc[(rows['z'] > 0)]
    rows4 = rows4.append(rows.loc[(rows['y'] == -102.5) & (rows['z'] == 0)])
    rows4 = rows4.append(rows.loc[(rows['y'] == 102.5) & (rows['z'] == 0)])
    rows4 = rows4.sort_values(by='y',ascending=False)
    avg4 = rows4.groupby(['y', 'z']).mean().reset_index()
    points4 = np.array([avg4['z'].to_numpy(), avg4['y'].to_numpy()]).T.reshape(-1, 1, 2)
    segments4 = np.concatenate([points4[:-1], points4[1:]], axis=1)
    lc4 = LineCollection(segments4, cmap='jet', norm=norm)
    lc4.set_array(avg4['avgMises'])
    lc4.set_linewidth(3)
    line = ax.add_collection(lc4)

    plt.grid()
    plt.axis('equal')
    ax.autoscale()
    plt.colorbar(line)
    #multiline(avg1['z'],avg1['y'],avg1['c'],cmap='jet')
    #plt.show()
    #newrows = rows.groupby(['y', 'z']).mean().reset_index()
    #plt.scatter(newrows['z'],newrows['y'],c=newrows['c'],cmap="jet")
    #plt.colorbar()
    ##minAvgMises = min(rows['avgMises'])
    plt.show()

    return

def twistAlongX(case):
    if case == "bending":
        Df = pd.read_csv('validation_processed/bending_processed.csv', header=0,index_col=0)
    elif case == "jam_bent":
        Df = pd.read_csv('validation_processed/jam_bent_processed.csv', header=0,index_col=0)
    elif case == "jam_straight":
        Df = pd.read_csv('validation_processed/jam_straight_processed.csv', header=0,index_col=0)
    else:
        print("Case not found")
        return
    yHinge = 0
    zHinge = 0
    yLE = 0
    zLE = 102.5

    rowsHinge = Df.loc[(Df['z'] == zHinge) & (Df['y'] == yHinge)]
    u2HingeAvg = rowsHinge.groupby(['x']).mean().reset_index()
    #u2HingeAvg = rowsHinge.drop_duplicates(subset=['x'])
    u2HingeAvg = u2HingeAvg.sort_values(by=['x'])

    rowsLE = Df.loc[(Df['z'] == zLE) & (Df['y'] == yLE)]
    u2LEAvg = rowsLE.groupby(['x']).mean().reset_index()
    #u2LEAvg = rowsLE.drop_duplicates(subset=['x'])
    u2LEAvg = u2LEAvg.sort_values(by=['x'])


    difY = u2HingeAvg['u2Loc1'].to_numpy()-u2LEAvg['u2Loc1'].to_numpy()
    difZ = zLE
    twist = np.arctan2(difY,difZ)

    fig, ax = plt.subplots()
    plt.suptitle("Twist angle VS distance in x-direction for case: %s" % (case), fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    ax.set_ylabel('Theta [rad]', fontsize=12.0)
    ax.set_xlabel('x [mm]', fontsize=12.0)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    plt.plot(u2HingeAvg['x'], twist)
    plt.grid()
    ax.autoscale()
    plt.show()

    return


#deflectionAlongX(100,0,"bending")
#maxStressAlongX("bending")
#stressCrossSection(1074,"bending")
twistAlongX("jam_bent")
