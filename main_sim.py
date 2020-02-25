# -*- coding: utf-8 -*-
"""
SVV 2020- Structural Analysis Assignment

Simulation for stress and deflection in A320 Airleron


@author: vladg
@version: 24-02-#3

"""
import numpy as np
import scipy as sp
from integration import *
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from calcNormStress import *
from EqSolvV2 import *
from verification import *


plt.close('all')
np.set_printoptions(precision=7)

#++++++++++++++++Global variables ++++++++++++++++


#++++++++++++++++++++++++++ Constants +++++++++++++++++++++++++++++++++++++++++++++++++++
g = 9.81 #m/s2


#+++++++++++++++++++++++++++ Main Simulation Functions ++++++++++++++++++++++++++++++++
def drawSection(ha,ca,stringer_posz,stringer_posy,Zcg,Zsc): #Verified by Vlad & Alberto!
    """Plots the cross-section."""
    fig,ax = plt.subplots()
    plt.title("Cross section")

    an = np.linspace(np.pi/2, np.pi*3/2, 100)
    plt.plot(ha/2+ha/2 * np.cos(an), ha/2 * np.sin(an))
    lines = [[(ha/2, ha/2), (ha/2, -ha/2)],
             [(ha/2,ha/2),(ca,0)],[(ha/2,-ha/2),(ca,0)]]
    lc = mc.LineCollection(lines, linewidths=2)
    ax.add_collection(lc)
    #Draw the centroid and Shear Centre
    plt.plot(-Zcg, 0, "og",label='Centroid')
    plt.plot(-Zsc,0,"or",label="S.C.")
    plt.grid()
    plt.axis('equal')
    plt.scatter(stringer_posz,stringer_posy,label='Stiffener')
    plt.legend()
    ax.autoscale()

#++++++++++++++++++++ Stringer Position (from TE to LE and back in c.c. order)++++++++++
def calcStPose(ha,ca,nst):   #Verified by Vlad!
    """Calculates the (y,z) position of the stringers.
    Input: height aileron, chord, number stringers
    Output: 2xnst array with the (y,z) coordinates"""

    if(nst <=0 or nst%2==0):
        print("Invalid input entered!")
        return 0

    pos = np.zeros((2,nst))
    #pos[:,int(nst/2)+1] = [0,0]
    semiCircum = calcCircum(ha,ca)
    unit     = semiCircum/ (int(nst))
    stCircle = int((ha/2 * np.pi/2 / unit))

    #Going around the semic cirle
    for i in range(1,stCircle+1):
        pos[:,int(nst/2)+1 - i - 1] = [ha/2 * np.sin(i*unit/(ha/2)),-ha/2*(1 - np.cos(i*unit/(ha/2)))]
        pos[:,int(nst/2)+1 + i -1 ] = [-ha/2 *np.sin(i*unit/(ha/2)),-ha/2*(1 - np.cos(i*unit/(ha/2)))]
     #Going along the skin
    if (2*stCircle +1 < nst):
        left   = (stCircle+1)*unit - ha/2 * np.pi/2     #left distance from the unit

        leftSt = int(nst/2) - stCircle     #nr of stringers left to palce outside the circle
        alfa   = np.arctan2(ha/2,ca-ha/2)   #the slope angle of the straight part

        leftTE = unit/2
        for i in range(leftSt):
             pos[:,i]  = [(leftTE+unit*i)*np.sin(alfa),-ca + (leftTE+unit*i)*np.cos(alfa)]
             pos[:,len(pos[1,:])-i-1] = [-(leftTE+unit*i)*np.sin(alfa),-ca + (leftTE+unit*i)*np.cos(alfa)]
        return pos
    else:
        return pos[1]


#+++++++++++++++++++++++++++++++ Geomtry parameters +++++++++++++++++++++++++++++++++++++++++
def calcCircum(ha,ca):
    #Stick to the name given in the flow chart OR
    # cleary write what you cahnged
    """Input: ha, ca
       Output: lCirc = length of the section circumference
    """
    return np.pi*ha/2 + 2*(np.sqrt((ca-ha/2)**2 + (ha/2)**2))

def calcCellArea(ha,ca):
    """Calculate area of left and right cells respectively"""
    # Input: ha, ca
    # Output: A1, A2
    A1 = (np.pi/2)*(ha/2)**2
    A2 = (ca-ha/2)*(ha/2)
    return A1,A2

def calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst):  #Verified by Vlad!
    """Return the z coordinate of the centroid."""
    stArea = calcStArea(tst,hst,wst)
    stPos  = calcStPose(ha,ca,nst)

    plateYLength = ha / 2
    plateZLength = ca - ha / 2

    sumStAreaZ   = np.sum(stPos[1,:]*stArea)
    sumAreaZ     = np.pi*tsk*(ha/2) * (-ha/2+(2/np.pi)*(ha/2)) + np.linalg.norm([plateYLength,plateZLength])*tsk*2 * (-ha/2 - plateZLength/2)+ ha*tsp*(-ha/2) + sumStAreaZ

    sumArea      = np.pi*tsk*(ha/2) + np.linalg.norm([plateYLength,plateZLength])*tsk*2 + ha*tsp + stArea*nst

    zCentroid    = sumAreaZ/sumArea

    return zCentroid        #note Zcg is negative

def calcInertia(Ca,H,Tsk,Tsp,Tst,Ast,Zcg,StPos):
    """Calculates the moment of inertia for Izz and Iyy and outputs in this order"""
    """Method is verified and fully correct"""


    Lsk = np.linalg.norm([H/2,Ca-H/2])     #Length of the slanted skin


    # ------------------------   Izz   ------------------------
    # I_zz consists of 4 parts: spar (1), skin plates (2), skin semicircular (3), stiffeners (4)

    I_zz1 = 1/12*Tsp*H**3
    Beta_plate = np.arctan2((H/2),(Ca-H/2))
    I_zz2 = 1/12*Tsk*(2*Lsk)**3*(np.sin(Beta_plate))**2
    I_zz3 = 2*(np.pi/4) *Tsk * ((H/2)**3) #0.5*1/8*np.pi*Tsk*H**3
    I_zz4 = Ast*sum(StPos[0,:]**2) #calcStPos gives list of coordinates (y,z)

    Izz = I_zz1 + I_zz2+I_zz3+I_zz4
    #print(I_zz1,I_zz2,I_zz3,I_zz4)

    # ------------------------   Iyy   -------------------------
    # I_yy consists of 4 parts: skin plates (1), skin semicircular (2), stiffeners (3), spar (4)
    # the MoI of the thinwalled semicircle about diameter is r^3*t*pi/4 (calculated by hand,verified)
    # the MoI of the thinwalled semicircle about cg is r^3*t*(pi/2 - 4/pi) (calculated by hand)

    I_yy_plate = 1/12*Tsk*(Lsk)**3*(np.cos(Beta_plate))**2 + Lsk*Tsk*(-H/2-0.5*(Ca-H/2)-Zcg)**2 #note the plus before Zcg because Zcg is negative itself
    I_yy1 = 2*I_yy_plate
    I_yy2 = (H/2)**3*Tsk*(np.pi/2 - 4/np.pi) + np.pi*H/2*Tsk * ((-H/2+H/np.pi)-Zcg)**2
    I_yy3 = Ast*sum((StPos[1,:]-Zcg)**2)   #calcStPos gives (y,z)
    I_yy4 = H*Tsp*(-H/2-Zcg)**2  #only steiner term due to thin walled approx (difference 5*e-10)

    Iyy = I_yy1+I_yy2+I_yy3+I_yy4
    #print(I_yy1,I_yy2,I_yy3,I_yy4)

    return Izz, Iyy

def calcTorsionStiffness(ha,ca,tsk,tsp,G):
    # Calculate torsional stiffness
    T = 1
    q1,q2,q3,q4,q01,q02,dthetadx = calcShFlowTorque(ha,ca,tsk,tsp,G,T)
    J = T/(G*dthetadx)
    return J

def calcStArea(tst, hst, wst):   #Verified by Vlad!
    #Calculates area of stringer in m^2
    StArea = tst * (hst + wst)
    return StArea

#+++++++++++++++++++++++++ Stress calculations +++++++++++++++++++++++++++++++++++++++++
def calcShCenter(ha,ca,tsk,tsp, tst, hst, wst,nst,n1,n2,n3,n4):
    """Calculates Z coordinate of shear center (according to the used coordiante system)"""

    #---------------------------------------------------
    # Calculates Z coordinate of shear center
    #
    # output:   description:                                                                                    type:
    # zShear    z-coordinate of the shear center, origin on leading edge, trailing edge is negative direction   value
    #
    #---------------------------------------------------

    Sy = 1
    Sz = 0
    # open section shear flow
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4 = calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy,n1,n2,n3,n4)

    # dimensions of upper and lower plate
    plateYLength = ha / 2
    plateZLength = ca - ha / 2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))

    #Let's draw some stuff...
    #drawGraph(sVec1,q1)

    # calculate moment around center of semi circle to find moment arm
    zeta = integrationArray(q4, sVec4[0],sVec4[-1],n4)*(ha / 2) + integrationArray(q1 ,sVec1[0],sVec1[-1],n1) * (plateZLength / plateLength) * (ha / 2) + integrationArray(q3, sVec3[0],sVec3[-1],n3) * (plateZLength / plateLength) * (ha / 2)

    # old zeta using scipy.integrate.simps()
    # zeta = sp.integrate.simps(q4, sVec4)*(ha / 2) + sp.integrate.simps(q1, sVec1) * (plateZLength / plateLength) * (ha / 2) + sp.integrate.simps(q3, sVec3) * (plateZLength / plateLength) * (ha / 2)

    zShear = -zeta - (ha / 2)

    return zShear


def calcShFlowTorque(ha,ca,tsk,tsp,G,T):
    #---------------------------------------------------
    # Calculate shear flow in both cells caused by torque
    #
    # output:   description:                                                                            type:
    # q1        shear flow in upper plate, positive in direction of increasing y and increasing z       value
    # q2        shear flow in spar, positive in direction of decreasing y                               value
    # q3        shear flow in lower plate, positive in direction of increasing y and decreasing z       value
    # q4        shear flow in semicircle, positive in direction of decreasing y (c.c)                        value
    #
    # q1, q2, q3, q4 are singular values that need to be added to the shear flow caused by shear forces acting through the shear center (calcShFlow)
    #
    #---------------------------------------------------

    # geometry calculations
    A1,A2 = calcCellArea(ha,ca)
    plateYLength = ha / 2
    plateZLength = ca - ha / 2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))

    # setting up matrix equations
    a11 = (1/(2*A1*G)) * ((np.pi*ha)/(2*tsk)+(ha/tsp))
    a12 = -(1/(2*A1*G)) * (ha/tsp)
    a13 = -1
    a21 = -(1/(2*A2*G)) * (ha/tsp)
    a22 = (1/(2*A2*G)) * ((2*plateLength)/(tsk)+(ha/tsp))
    a23 = -1
    a31 = 2*A1
    a32 = 2*A2

    a33 = 0  #only two closed sections, so this must twist rate
    A = np.matrix([[a11,a12,a13],
                   [a21,a22,a23],
                   [a31,a32,a33]])

    b = np.array([0,
                  0,
                  T])

    # x[0]: q1 (constant shear flow left cell)
    # x[1]: q2 (constant shear flow right cell)
    # x[2]: dTheta/dx (twist rate)
    # q1 and q2 are positive counter clockwise
    x = np.linalg.solve(A,b)

    q1 = x[1]
    q2 = x[1]-x[0]
    q3 = x[1]
    q4 = x[0]

    return q1,q2,q3,q4,x[0],x[1],x[2]


def calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy,n1,n2,n3,n4):
    print("Calculating shear flows...")

    zCentroid   = calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst)
    stArea      = calcStArea(tst,hst,wst)

    #---------------------------------------------------
    # Calculate shear flow, under condition that the shear forces act through the shear center
    #
    # output:   description:                                                                            type:
    # q1        shear flow in upper plate, positive in direction of increasing y and increasing z       array[n1+1 , 1]
    # q2        shear flow in spar, positive in direction of decreasing y                               array[n2+1 , 1]
    # q3        shear flow in lower plate, positive in direction of decreasing y and decreasing z       array[n3+1 , 1]
    # q4        shear flow in semicircle, positive in direction of decreasing y                         array[n4+1 , 1]
    #
    # q1, q2, q3, q4 are singular values that need to be added to the shear flow caused by shear forces acting through the shear center (calcShFlow)
    #
    #---------------------------------------------------

    # geometry calculations
    zCentroid = calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst)
    stArea = calcStArea(tst,hst,wst)

    stringerPos = calcStPose(ha, ca, nst)
    Izz,Iyy     = calcInertia(ca, ha, tsk, tsp, tst, stArea, zCentroid, stringerPos)

    # calculate stringer coordinates w.r.t. centroid
    stringerPosCentroid = stringerPos
    stringerPosCentroid[1, :] = stringerPosCentroid[1, :] - zCentroid

    # upper and lower flat plate lengths
    plateYLength = ha/2
    plateZLength = ca-ha/2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))

    #---------------------------------------------------
    # (1) top flat plate : starting from (y = 0, z = -ca)
    #---------------------------------------------------

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

    #---------------------------------------------------
    # (2) spar plate : continuing from (1)
    #---------------------------------------------------
    sVec2 = np.linspace(0,ha,n2+1)
    yVec2 = ha/2-sVec2
    zVecCentroid2 = np.ones(len(sVec2))*(-zCentroid -ha/2)
    ds2 = ha/n2
    intYVec2,intZVec2 = calcIntegralArray(zVecCentroid2,yVec2,sVec2,n2,tsp)
    qs2 = -(Sz/Iyy)*intZVec2 - (Sy/Izz)*intYVec2
    # add last value from qs1
    qs2 += qs1[-1]

    #---------------------------------------------------
    # (4) semi-circular area : starting from (y = ha/2, z = -ha/2)
    #---------------------------------------------------

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


    #---------------------------------------------------
    # (3) lower flat plate : continuing from (2) and (4)
    #---------------------------------------------------

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
    # add last value from qs2 to qs3
    qs3 += qs2[-1]
    # add last value from qs4 to qs3
    qs3 += qs4[-1]

    #---------------------------------------------------
    # calculating constant shear flow
    #---------------------------------------------------

    # old, using scipy.integrate.simps()
    # rhs1 = -(-sp.integrate.simps(qs2, sVec2) / tsp + sp.integrate.simps(qs4, sVec4) / tsk)
    # rhs2 = -(sp.integrate.simps(qs1, sVec1) / tsk + sp.integrate.simps(qs2, sVec2) / tsp + sp.integrate.simps(qs3,sVec3) / tsk)

    # setting up matrix equations
    rhs1 = -(-integrationArray(qs2, sVec2[0],sVec2[-1],n2) / tsp + integrationArray(qs4, sVec4[0],sVec4[-1],n4) / tsk)
    rhs2 = -(integrationArray(qs1, sVec1[0],sVec1[-1],n1) / tsk + integrationArray(qs2, sVec2[0],sVec2[-1],n2) / tsp + integrationArray(qs3, sVec3[0],sVec3[-1],n3) / tsk)
    A = np.matrix([[(np.pi * ha) / (2 * tsk) + (ha / tsp), -(ha / tsp)],
                   [-(ha / tsp), (2 * plateLength) / (tsk) + (ha / tsp)]])
    b = np.array([rhs1, rhs2])

    # x[0]: q1 (constant shear flow of left cell)
    # x[1]: q2 (constant shear flow of right cell)
    # q1 and q2 are positive counter-clockwise
    x = np.linalg.solve(A, b)

    # adding constant shear flow to open section shear flow to get final shear flow
    q1 = qs1 + x[1]
    q2 = qs2 + x[1] - x[0]
    q3 = qs3 + x[1]
    q4 = qs4 + x[0]


    q1 = np.vstack((yVec1,zVecCentroid1,q1))
    q2 = np.vstack((yVec2,zVecCentroid2,q2))
    q3 = np.vstack((yVec3,zVecCentroid3,q3))
    q4 = np.vstack((yVec4,zVecCentroid4,q4))


    return q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4

#++++++++++++++++++++++++++++ Draw shear flow +++++++++++++++++++++++++++++++++++++++
def drawGraph(x,y):

    """Helper function to draw graphs"""

    print("Draw shear flow distribution...")

    #---------------------------------------------------
    # Helper function to draw graphs
    #
    #---------------------------------------------------

    fig, ax = plt.subplots()
    plt.title("Shear Flow")
    plt.plot(x,y)
    plt.grid(True)
    ax.autoscale()
    plt.show()


def calcIntegralArray(z,y,s, n,t):
    #---------------------------------------------------
    # Evaluate integrals needed in calcShFlow()
    #
    # output:       description:                                type:
    # t*intYVec     integral of (t*y*ds) from s[0] to s[-1]     array[n+1 , 1]
    # t*intZVec     integral of (t*z*ds) from s[0] to s[-1]     array[n+1 , 1]
    #
    # Assuming t is constant in the section that is integrated from s[0] to s[-1]
    #
    #---------------------------------------------------

    intYVec = np.array([0])
    intZVec = np.array([0])
    for i in range(1, n+1):
        # old, using sp.integrate.simps()
        # integralY = sp.integrate.simps(y[:i + 1], s[:i + 1])
        integralY = integrationArray(y[:i+1],s[0],s[i],i)
        intYVec = np.append(intYVec, integralY)

        # old, using sp.integrate.simps()
        # integralZ = sp.integrate.simps(z[:i + 1], s[:i + 1])
        integralZ = integrationArray(z[:i+1],s[0],s[i],i)
        intZVec = np.append(intZVec, integralZ)

    return t*intYVec,t*intZVec


def addStringerContribution(integrated,yVec,zVec,ds,n,stringerPos,stArea,direction):
    #---------------------------------------------------
    # Add stringer contribution to shear flow
    #
    # output:           description:                                                                    type:
    # newIntegrated     integral of (t*y*ds) OR (t*z*ds) from s[0] to s[-1] + stringer contribution     array[n+1 , 1]
    #
    # Stringers are booms with point area B, added contribution is B*y OR B*z
    #
    #---------------------------------------------------
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


def calcCellArea(ha,ca):
    """Calculate area of left and right cells respectively"""
    # Input: ha, ca
    # Output: A1, A2
    A1 = (np.pi/2)*(ha/2)**2
    A2 = (ca-ha/2)*(ha/2)
    return A1,A2


def calcDist(y1,z1,y2,z2):
    #---------------------------------------------------
    # Calculate distance between two points (Euclidean norm)
    #
    # output:       description:                                type:
    # dist          distance between the two given points       value
    #
    #---------------------------------------------------
    dist = np.sqrt(np.power(y1-y2, 2) + np.power(z1-z2, 2))

    return dist



def calcNormStress(y,z,Iyy,Izz,MomentArray):
    #@author: alberto
    # Only 2nd and 3rd row of MomentArray are used


    # Append spar cap locations to StPos
    #Moment Array =[[T][My][Mz]]

    # ---------------- Sigma -------------------------
    # the goal is to create an array with only the maximum stress and its location for each section
    # This array is MaxStresses

    # print(MomentArray)
    sigma_y =  float(MomentArray[1,:])*z/Iyy     # compute stress
    sigma_z = -float(MomentArray[2,:])*y/Izz
    sigma   = sigma_y + sigma_z
    return sigma

def calcTau(T,Sy,Sz,ha,ca,tsk,tsp, tst, hst, wst,nst,G,n1,n2,n3,n4):
    """Returns tau at the the location x"""
    q1,q2,q3,q4,s1,s2,s3,s4  = calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy,n1,n2,n3,n4)
    q1t,q2t,q3t,q4t,x1,x2,x3 = calcShFlowTorque(ha,ca,tsk,tsp,G,T)

    q1[2,:]+= q1t
    q2[2,:]+= q2t
    q3[2,:]+= q3t
    q4[2,:]+= q4t

    tau1 = q1
    tau1[2,:] = q1[2,:]/tsk

    tau2 =q2
    tau2[2,:] = q2[2,:]/tsp    #spar top to bottom
    tau3=q3
    tau3[2,:] = q3[2,:]/tsk
    tau4=q4
    tau4[2,:] = q4[2,:]/tsk

    tau = np.hstack((tau1,tau4,tau3,tau2))   #cc around the section and than the spa
    # print(np.shape(tau))
    return tau

def VonMisses(sigma,tau):
    """Returns the combined loading stress according to von Misses yeild criterion.
    Input:
           sigma = [[z],[y],[sigma_x]]
           tau   = ask Danny
            """
    A = np.power(sigma,2)
    B = np.power(tau,2)

    #sigm zz and sigma_yy are zero

    return np.sqrt(A/2 + 3*B)



#++++++++++++++++++++++++++ Shear & Moment ++++++++++++++++++++++++++++++++++++++++++++++++
def genMoments(x):
    """Imports the T,Mx,Mz function from the equation solver and calculatates all the moment
    as a function of x"""
    Tval =  T(x)
    Myval = My(x)
    Mzval = Mz(x)

    return np.array([[Tval],[Myval],Mzval])

def genVM(x,craft):
    """Integrates toghether the moment evaluation, sigma is taken.
    """
    MomentArray = genMoments(x) #[[T],[My],[Mz]]

    #sigma       = calcNormStress(Ha, ca, Zcg, Iyy, Izz, MomentArray, npoints)

    tau   = calcTau(float(MomentArray[0]),0,0,craft.ha,craft.ca,craft.tsk,craft.tsp, craft.tst, craft.hst, craft.wst,craft.nst,craft.G,craft.n1,craft.n2,craft.n3,craft.n4)#Change this!!

    y=    tau[0,:]
    z=    tau[1,:]
    tau = tau[2,:]


    sigma = calcNormStress(y,z,craft.Izz,craft.Iyy,MomentArray)
    print("Sigma:\n",sigma)
    VM    = VonMisses(sigma,tau)

    return sigma,tau,VM,y,z


#++++++++++++++++++++++++++++++ Aircraft Class ++++++++++++++++++++++++++++++++++++++++++++
class Aircraft:
    def __init__(self,name):
        if name=="A320" or name=="a320":
            self.n1 = 20
            self.n2 = 20
            self.n3 = 20
            self.n4 = 20
            self.la  = 2.771          #m
            self.ca  = 0.547
            self.ha  = 0.225
            self.tsk = 0.0011
            self.tst = 0.0012
            self.wst = 0.020
            self.hst = 0.015
            self.nst = 17
            self.tsp = 0.0029
            self.Ast = calcStArea(self.tst,self.hst,self.wst)
            self.theta = np.radians(26)  #rad
            self.E     = 73.1*10**9     #aluminium 2024-T3
            self.G     = 28*10**9       #aluminium 2024-T3
            self.Zcg   = calcCentroid(self.ha,self.ca,self.tsk,self.tsp,self.tst,self.hst,self.wst,self.nst)
            self.StPos = calcStPose(self.ha, self.ca, self.nst)
            self.Iyy,self.Izz = calcInertia(self.ca,self.ha,self.tsk,self.tsp,self.tst,self.Ast,self.Zcg,self.StPos)

            
#++++++++++++++++++++++++++++ Main +++++++++++++++++++++++++++++++++++++++++++++++++++

def main():
    craft = Aircraft("A320")

    # print("Circumference: \n",calcCircum(craft.ha,craft.ca))

    # stArea = calcStArea(craft.tst,craft.hst,craft.wst)
    # #print("Stringer Area is:\n",stArea)

    # A1,A2 = calcCellArea(craft.ha,craft.ca)
    # #print("Cell areas are:\n",A1,A2)

    # pos = calcStPose(craft.ha,craft.ca,craft.nst)
    # #print("Stringers (y,z) are:\n",pos)

    Zcg = calcCentroid(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst)
    # #print("Centroid z-coordinate is:\n", Zcg)

    # Izz,Iyy = calcInertia(craft.ca,craft.ha,craft.tsk,craft.tsp,craft.tst,craft.Ast,Zcg,pos)
    # #print("Izz and Iyy:\n",Izz, Iyy)

    # q = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,1,0, discret.n1,discret.n2,discret.n3,discret.n4)
    # #print("Shear flows are:\n", q)

    # Zsc = calcShCenter(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,discret.n1,discret.n2,discret.n3,discret.n4)
    # #print("Shear center z-coordinate is:\n", zShear)

    # vm  = VonMisses(np.array([[0],[0],[0]]),q)
    # #print("Von Misses stress are:\n",vm)

    #drawSection(craft.ha,craft.ca,-pos[1,:],-pos[0,:],Zcg,Zsc)

    # J = calcTorsionStiffness(craft.ha, craft.ca, craft.tsk, craft.tsp, craft.G)
    # #print("J:\n",J)

    # zShear = calcShCenter(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,discret.n1,discret.n2,discret.n3,discret.n4)
    # #print("Shear center z-coordinate is:\n", zShear)

    #Pick a spanwise x location
    x1 =0.5
    n = craft.n1 + craft.n2 + craft.n3 +craft.n4  #keep it odd
    sigma,tau,VM,y,z = genVM(x1,craft)

    plt.scatter(-Zcg-z,y,c=VM, vmin=min(VM), vmax=max(VM), s=20, cmap="jet")
    print(VM/(10**6))

    #drawSection(craft.ha,craft.ca,-pos[1,:],-pos[0,:],Zcg,Zsc)

    #shearFlowGraph2(craft)



if __name__ == "__main__":
    main()
    plt.show()
