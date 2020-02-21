"""
Created on 20-02 16:50
Shear flow and shear center
@author: dannywhuang
@version: 21-02
"""

import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from main_sim import *
from integration import *

def calcShCenter(ha,ca,tsk,tsp, tst, hst, wst,nst,n1,n2,n3,n4):
    """Calculates Z coordinate of shear center"""
    Sy = 1
    Sz = 0
    # open section shear flow
    qs1,qs2,qs3,qs4,sVec1,sVec2,sVec3,sVec4 = calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy,n1,n2,n3,n4)
    # dimensions of upper and lower plate
    plateYLength = ha / 2
    plateZLength = ca - ha / 2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))
    # calculate const. sh flow
    #rhs1 = -(-sp.integrate.simps(qs2, sVec2) / tsp + sp.integrate.simps(qs4, sVec4) / tsk)
    #rhs2 = -(sp.integrate.simps(qs1, sVec1) / tsk + sp.integrate.simps(qs2, sVec2) / tsp + sp.integrate.simps(qs3,sVec3) / tsk)
    rhs1 = -(-integrationArray(qs2, sVec2[0],sVec2[-1],n2) / tsp + integrationArray(qs4, sVec4[0],sVec4[-1],n4) / tsk)
    rhs2 = -(integrationArray(qs1, sVec1[0],sVec1[-1],n1) / tsk + integrationArray(qs2, sVec2[0],sVec2[-1],n2) / tsp + integrationArray(qs3, sVec3[0],sVec3[-1],n3) / tsk)

    A = np.matrix([[(np.pi * ha) / (2 * tsk) + (ha / tsp), -(ha / tsp)],
                   [-(ha / tsp), (2 * plateLength) / (tsk) + (ha / tsp)]])
    b = np.array([rhs1, rhs2])
    # x[0] is q,0 of left cell, x[1] is q,0 of right cell
    x = np.linalg.solve(A, b)
    # add const. shear flow to open section shear flow to get final shear flow
    q1 = qs1 + x[1]
    q2 = qs2 + x[1] - x[0]
    q3 = qs3 + x[1]
    q4 = qs4 + x[0]

    drawGraph(sVec1,q1)
    # calculate moment around center of semi circle to find moment arm
    #zeta = sp.integrate.simps(q4, sVec4)*(ha / 2) + sp.integrate.simps(q1, sVec1) * (plateZLength / plateLength) * (ha / 2) + sp.integrate.simps(q3, sVec3) * (plateZLength / plateLength) * (ha / 2)
    zeta = integrationArray(q4, sVec4[0],sVec4[-1],n4)*(ha / 2) + integrationArray(q1 ,sVec1[0],sVec1[-1],n1) * (plateZLength / plateLength) * (ha / 2) + integrationArray(q3, sVec3[0],sVec3[-1],n3) * (plateZLength / plateLength) * (ha / 2)
    zShear = -zeta - (ha / 2)
    return zShear


def calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy,n1,n2,n3,n4):
    """Calculate open section shear flow"""

    zCentroid = calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst)
    stArea = calcStArea(tst,hst,wst)
    stringerPos = calcStPose(ha, ca, nst)
    Izz,Iyy = calcInertia(ca, ha, tsk, tsp, tst, stArea, zCentroid, stringerPos)

    stringerPosCentroid = stringerPos
    # calculate stringer coordinates w.r.t. centroid
    stringerPosCentroid[1, :] = stringerPosCentroid[1, :] - zCentroid

    # upper and lower flat plate lengths
    plateYLength = ha/2
    plateZLength = ca-ha/2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))

    ###### (1) top flat plate : starting from (y = 0, z = -ca)
    # filter stringers of top flate plate
    stringerPosFilt1 = stringerPosCentroid[:, (stringerPosCentroid[0, :] >= 0) & (stringerPosCentroid[1, :] <= -(ha / 2)-zCentroid)]
    # sVec is the vector that contains s ----- yVec is the vector that contains y ----- zVecCentroid is the vector that contains z w.r.t. centroid
    sVec1 = np.linspace(0,plateLength,n1+1)
    yVec1 = (plateYLength/plateLength) * sVec1
    zVecCentroid1 = -zCentroid -ca + (plateZLength/plateLength) * sVec1
    # ds1 is the grid spacing
    ds1 = plateLength/n1
    # intYVec and intZVec are the integrals w.r.t. y and z respectively
    intYVec1,intZVec1 = calcIntegralArray(zVecCentroid1,yVec1,sVec1,n1,tsk)
    # stringer contribution is added afterwards
    intYVec1 = addStringerContribution(intYVec1,yVec1,zVecCentroid1,ds1,n1,stringerPosFilt1,stArea,"Y")
    intZVec1 = addStringerContribution(intZVec1,yVec1,zVecCentroid1,ds1,n1,stringerPosFilt1,stArea,"Z")
    qs1 = -(Sz/Iyy)*intZVec1 - (Sy/Izz)*intYVec1

    ###### (2) spar plate : continuing from (1)
    sVec2 = np.linspace(0,ha,n2+1)
    yVec2 = ha/2-sVec2
    zVecCentroid2 = np.ones(len(sVec2))*(-zCentroid -ha/2)
    ds2 = ha/n2
    intYVec2,intZVec2 = calcIntegralArray(zVecCentroid2,yVec2,sVec2,n2,tsp)
    qs2 = -(Sz/Iyy)*intZVec2 - (Sy/Izz)*intYVec2
    # add last value from qs1
    qs2 += qs1[-1]

    ###### (4) semi-circular area : starting from (y = ha/2, z = -ha/2)
    # filter stringers of semicircle
    stringerPosFilt4 = stringerPosCentroid[:, (stringerPosCentroid[1, :] >= (-ha / 2)-zCentroid)]
    thetaVec4 = np.linspace(0, np.pi, n4 + 1)
    sVec4 = thetaVec4*(ha/2)
    yVec4 =(ha / 2) * np.cos(thetaVec4)
    zVecCentroid4 =(-zCentroid - (ha/2) + (ha/2)*np.sin(thetaVec4))
    # yVecR is yVec multiplied by the radius ----- zVecCentroidR is zVecCentroid multiplied by the radius
    yVecR4 = yVec4 *  (ha / 2)
    zVecCentroidR4 = zVecCentroid4 * (ha/2)
    ds4 = (np.pi*(ha/2)) / n4
    intYVec4, intZVec4 = calcIntegralArray(zVecCentroidR4, yVecR4, thetaVec4, n4, tsk)
    intYVec4 = addStringerContribution(intYVec4, yVec4, zVecCentroid4, ds4, n4, stringerPosFilt4, stArea, "Y")
    intZVec4 = addStringerContribution(intZVec4, yVec4, zVecCentroid4, ds4, n4, stringerPosFilt4, stArea, "Z")
    qs4 = -(Sz / Iyy) * intZVec4 - (Sy / Izz) * intYVec4

    ##### (3) lower flat plate : continuing from (2) and (4)
    # filter stringers of lower flate plate
    stringerPosFilt3 = stringerPosCentroid[:,(stringerPosCentroid[0, :] <= 0) & (stringerPosCentroid[1, :] <= (-ha / 2) - zCentroid)]
    sVec3 = np.linspace(0, plateLength, n3 + 1)
    yVec3 = (-ha / 2) + (plateYLength / plateLength) * sVec3
    zVecCentroid3 = (-zCentroid - (ha / 2)) - (plateZLength / plateLength) * sVec3
    ds3 = plateLength / n3
    intYVec3, intZVec3 = calcIntegralArray(zVecCentroid3, yVec3, sVec3, n3, tsk)
    intYVec3 = addStringerContribution(intYVec3, yVec3, zVecCentroid3, ds3, n3, stringerPosFilt3, stArea, "Y")
    intZVec3 = addStringerContribution(intZVec3, yVec3, zVecCentroid3, ds3, n3, stringerPosFilt3, stArea, "Z")
    qs3 = -(Sz / Iyy) * intZVec3 - (Sy / Izz) * intYVec3

    #add last value from qs2 to qs3
    qs3 += qs2[-1]
    #add last value from qs4 to qs3
    qs3 += qs4[-1]

    return qs1,qs2,qs3,qs4,sVec1,sVec2,sVec3,sVec4


def drawGraph(x,y):
    """Helper function to draw graphs"""
    fig, ax = plt.subplots()
    plt.title("Shear Flow")
    plt.plot(x,y)
    plt.grid()
    ax.autoscale()
    plt.show()


def calcIntegralArray(z,y,s, n,t):
    """Used in shear flow calculations"""
    intYVec = np.array([0])
    intZVec = np.array([0])
    for i in range(1, n+1):
        #integralY = sp.integrate.simps(y[:i + 1], s[:i + 1])
        integralY = integrationArray(y[:i+1],s[0],s[i],i)
        intYVec = np.append(intYVec, integralY)

        #integralZ = sp.integrate.simps(z[:i + 1], s[:i + 1])
        integralZ = integrationArray(z[:i+1],s[0],s[i],i)
        intZVec = np.append(intZVec, integralZ)

    return t*intYVec,t*intZVec


def addStringerContribution(integrated,yVec,zVec,ds,n,stringerPos,stArea,direction):
    """Used in shear flow calculations"""
    used = np.array([])
    newIntegrated = integrated
    for i in range(0,n+1):
        for j in range(0,stringerPos.shape[1]):
            dist = calcDist(stringerPos[0,j],stringerPos[1,j],yVec[i],zVec[i])
            if dist < ds and j not in used:
                used = np.append(used,j)
                if direction=="Y":
                    newIntegrated[i+1:] += stArea*stringerPos[0,j]
                elif direction == "Z":
                    newIntegrated[i+ 1:] += stArea * stringerPos[1, j]
    if len(used) != stringerPos.shape[1]:
        print("Warning addStringerContribution(): Not all stringers have been used")
        print("Amount of stringers not used:\n", stringerPos.shape[1] - used.size)

    return newIntegrated


def calcDist(y1,z1,y2,z2):
    """Euclidean norm"""
    return np.sqrt(np.power(y1-y2, 2) + np.power(z1-z2, 2))