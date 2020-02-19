# -*- coding: utf-8 -*-
"""
SVV 2020- Structural Analysis Assignment

Simulation for stress and deflection in A320 Airleron
19/02 morning session
@author: vladg
@version: 19-02-#1
"""
import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
plt.close('all')
np.set_printoptions(precision=3)

#++++++++++++++++++++++++++ Constants +++++++++++++++++++++++++++++++++++++++++++++++++++
g = 9.81 #m/s2


#+++++++++++++++++++++++++++ Main Simulation Functions ++++++++++++++++++++++++++++++++
def calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy):
    #open section shear flow

    Sz = 1
    Sy = 0
    #testcase
    #ha = 6
    #ca = 7
    #tsk = 1
    #tsp = 1
    #tst = 1
    #hst = 1
    #wst = 1
    #nst = 5
    #Sz = 1
    #Sy = 0
    #Izz = 1
    #Iyy = 1

    zCentroid = calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst)
    stArea = calcStArea(tst,hst,wst)
    stringerPos = calcStPose(ha, ca, nst)

    zCentroid = -0.21578
    Izz, Iyy = 1.28074*10**(-5),6.84137*10**(-5)
    #calcInertia(ca, ha, tsk, tsp, tst, stArea, zCentroid, stringerPos)
    print("zCentroid:\n",zCentroid)

    #stringerPos = np.array([[0,(ha/2)*np.sin(np.pi/4),-(ha/2)*np.sin(np.pi/4),1.5,-1.5],[0,-ha/2+(ha/2)*np.cos(np.pi/4),-ha/2+(ha/2)*np.cos(np.pi/4),-5,-5]])
    stringerPosCentroid = stringerPos
    stringerPosCentroid[1, :] = stringerPosCentroid[1, :] - zCentroid
    print("stringerPosCentroid:\n", stringerPosCentroid)

    plateYLength = ha/2
    plateZLength = ca-ha/2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))
    ### right cell
    ############ (1) top flat plate : starting from (y = 0, z = -ca)
    n1 = 1000
    #get stringers on top flate plate
    stringerPosFilt1 = stringerPosCentroid[:, (stringerPosCentroid[0, :] >= 0) & (stringerPosCentroid[1, :] <= -(ha / 2)-zCentroid)]
    sVec1 = np.linspace(0,plateLength,n1+1)
    yVec1 = (plateYLength/plateLength) * sVec1
    zVecCentroid1 = -zCentroid -ca + (plateZLength/plateLength) * sVec1
    ds1 = plateLength/n1
    intYVec1,intZVec1 = calcIntegralArray(zVecCentroid1,yVec1,sVec1,n1,tsk)
    intYVec1 = addStringerContribution(intYVec1,yVec1,zVecCentroid1,ds1,n1,stringerPosFilt1,stArea,"Y")
    intZVec1 = addStringerContribution(intZVec1,yVec1,zVecCentroid1,ds1,n1,stringerPosFilt1,stArea,"Z")
    qs1 = -(Sz/Iyy)*intZVec1 - (Sy/Izz)*intYVec1

    ############ (2) spar plate : continuing from (1)
    n2 = 1000
    #no stringers on spar duh
    sVec2 = np.linspace(0,ha,n2+1)
    yVec2 = ha/2-sVec2
    zVecCentroid2 = np.ones(len(sVec2))*(-zCentroid -ha/2)
    ds2 = ha/n2
    intYVec2,intZVec2 = calcIntegralArray(zVecCentroid2,yVec2,sVec2,n2,tsp)
    qs2 = -(Sz/Iyy)*intZVec2 - (Sy/Izz)*intYVec2
    #add last value from qs1
    qs2 += qs1[-1]
    
    ############ (3) lower flat plate : continuing from (2)
    n3 = 1000
    stringerPosFilt3 = stringerPosCentroid[:, (stringerPosCentroid[0, :] <= 0) & (stringerPosCentroid[1, :] <= (-ha / 2)-zCentroid)]
    sVec3 = np.linspace(0, plateLength, n3 + 1)
    yVec3 = (-ha/2) + (plateYLength / plateLength) * sVec3
    zVecCentroid3 = (-zCentroid-(ha/2)) - (plateZLength / plateLength) * sVec3
    ds3 = plateLength / n3
    intYVec3, intZVec3 = calcIntegralArray(zVecCentroid3, yVec3, sVec3, n3,tsk)
    intYVec3 = addStringerContribution(intYVec3, yVec3, zVecCentroid3, ds3, n3, stringerPosFilt3, stArea, "Y")
    intZVec3 = addStringerContribution(intZVec3, yVec3, zVecCentroid3, ds3, n3, stringerPosFilt3, stArea, "Z")
    qs3 = -(Sz / Iyy) * intZVec3 - (Sy / Izz) * intYVec3
    #add last value from qs2



    ###left cell
    ######## (4) semicircular arc starts here
    n4 = 1000
    stringerPosFilt4 = stringerPosCentroid[:, (stringerPosCentroid[1, :] >= (-ha / 2)-zCentroid)]
    thetaVec4 = np.linspace(0, np.pi, n4 + 1)
    sVec4 = thetaVec4*(ha/2)
    yVec4 =(ha / 2) * np.cos(thetaVec4)
    zVecCentroid4 =(-zCentroid - (ha/2) + (ha/2)*np.sin(thetaVec4))
    yVecR4 = yVec4 *  (ha / 2)
    zVecCentroidR4 = zVecCentroid4 * (ha/2)
    ds4 = (np.pi*(ha/2)) / n4
    intYVec4, intZVec4 = calcIntegralArray(zVecCentroidR4, yVecR4, thetaVec4, n4, tsk)
    intYVec4 = addStringerContribution(intYVec4, yVec4, zVecCentroid4, ds4, n4, stringerPosFilt4, stArea, "Y")
    intZVec4 = addStringerContribution(intZVec4, yVec4, zVecCentroid4, ds4, n4, stringerPosFilt4, stArea, "Z")
    qs4 = -(Sz / Iyy) * intZVec4 - (Sy / Izz) * intYVec4

    qs3 += qs2[-1]
    qs3 += qs4[-1]
    print("qs1", qs1)
    print("qs2", qs2)
    print("qs3", qs3)
    print("qs4", qs4)
    #calculate const. sh flow
    #integrate open section sh. flow for left cell
    rhs1 = -(-sp.integrate.simps(qs2,sVec2)/tsp+sp.integrate.simps(qs4,sVec4)/tsk)
    rhs2 = -(sp.integrate.simps(qs1,sVec1)/tsk+sp.integrate.simps(qs2,sVec2)/tsp+sp.integrate.simps(qs3,sVec3)/tsk)
    A = np.matrix([[(np.pi*ha)/(2*tsk)+(ha/tsp) ,   -(ha/tsp)],
                   [-(ha/tsp)                   ,   (2*plateLength)/(tsk)+(ha/tsp)]])
    b = np.array([rhs1,rhs2])
    #x[0] is q,0 of left cell, x[1] is q,0 of right cell
    x = np.linalg.solve(A,b)
    print(x)

    q1 = qs1 + x[1]
    q2 = qs2 + x[1] - x[0]
    q3 = qs3 + x[1]
    q4 = qs4 + x[0]
    print("q1",q1)
    print("q2",q2)
    print("q3",q3)
    print("q4",q4)
    return qs1,qs2,qs3,qs4



    #moment around hinge



def calcIntegralArray(z,y,s, n,t):
    intYVec = np.array([0])
    intZVec = np.array([0])
    for i in range(1, n+1):
        integralY = sp.integrate.simps(y[:i + 1], s[:i + 1])
        intYVec = np.append(intYVec, integralY)
        integralZ = sp.integrate.simps(z[:i + 1], s[:i + 1])
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
                print("yVec",yVec[i])
                print("zVec",zVec[i])
                if direction=="Y":
                    newIntegrated[j+1:] += stArea*stringerPos[0,j]
                elif direction == "Z":
                    newIntegrated[j + 1:] += stArea * stringerPos[1, j]
    print("used:\n",used)
    if len(used) != stringerPos.shape[1]:
        print("Warning addStringerContribution(): Not all stringers have been used")
    return newIntegrated

def calcDist(y1,z1,y2,z2):
    """Euclidean norm"""
    return np.sqrt(np.power(y1-y2, 2) + np.power(z1-z2, 2))


def drawSection(ha,ca,stringer_posz,stringer_posy,Zcg): #Verified by Vlad & Alberto!
    """Plots the cross-section."""
    fig,ax = plt.subplots()
    plt.title("Cross section")
    
    an = np.linspace(np.pi/2, np.pi*3/2, 100)
    plt.plot(ha/2+ha/2 * np.cos(an), ha/2 * np.sin(an))
    lines = [[(ha/2, ha/2), (ha/2, -ha/2)],
             [(ha/2,ha/2),(ca,0)],[(ha/2,-ha/2),(ca,0)]]
    lc = mc.LineCollection(lines, linewidths=2)
    ax.add_collection(lc)
    plt.plot(-Zcg, 0, "or",label='Centroid')
    
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
    semiCircum = calcCircum(ha,ca)/2

    unit     = semiCircum/ (int(nst/2) +1)

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
        lineSt = int((semiCircum - (ha/2 * np.pi/2) - left)/unit)

        leftTE = semiCircum- (stCircle+ lineSt)*unit
        for i in range(leftSt):
             pos[:,i]  = [(leftTE+unit*i)*np.sin(alfa),-ca + (leftTE+unit*i)*np.cos(alfa)]
             pos[:,len(pos[1,:])-i-1] = [-(leftTE+unit*i)*np.sin(alfa),-ca + (leftTE+unit*i)*np.cos(alfa)]
        return pos
    else:
        return pos[1]



def calcCircum(ha,ca):
    #Stick to the name given in the flow chart OR
    # cleary write what you cahnged
    """Input: ha, ca
       Output: lCirc = length of the section circumference
    """
    return np.pi*ha/2 + 2*(np.sqrt((ca-ha/2)**2 + (ha/2)**2))


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
    I_zz3 = 0.5*1/8*np.pi*Tsk*H**3
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


#++++++++++++++++++++++++++++++++ Numerical Integration +++++++++++++++++++++++++++++++++
def integration(function,n,a,b):

	zf = b

	deltaz = round(zf/n,8)
	zvector = np.array([0])
	h = deltaz/2
	primvector = np.array([0])
	valuez = 0

	for i in range(n):
		valuez += deltaz
		zvector = np.append(zvector,valuez)


	for j in range(n):
		value = h/3*(function(zvector[j]) + 4*function(zvector[j]+h) + function(zvector[j+1]))
		primvector = np.append(primvector,value)

	A = np.empty([n+1,n+1])
	for i in range(n+1):
		for j in range(n+1):
			A[i,j] = zvector[i]**j

	weights = np.dot(np.linalg.inv(A),np.transpose(primvector))
	total = np.sum(primvector)


	####  Verification ######
	#vervec = np.array([0])
	#verificationtotal = scp.integrate.simps(function(zvector),zvector)
	#for i in range(n):
		#vector = np.array([zvector[i],zvector[i+1]])
		#vervec = np.append(vervec,scp.integrate.simps(function(vector),vector))

	#plt.plot(zvector,vervec)
	#plt.plot(zvector,primvector)
	#plt.show()

	return weights,total



def calcStArea(Tst, Hst, Wst):   #Verified by Vlad!
  #Calculates area of stringer in m^2
    StArea = Tst * (Hst + Wst)
    return StArea

class Aircraft:
    def __init__(self,name):
        if name=="A320" or name=="a320":
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

#++++++++++++++++++++++++++++ Main +++++++++++++++++++++++++++++++++++++++++++++++++++

def main():

    craft = Aircraft("A320")
    print("Circumference: \n",calcCircum(craft.ha,craft.ca))

    stArea = calcStArea(craft.tst,craft.hst,craft.wst)
    print("Stringer Area is:\n",stArea)

    pos = calcStPose(craft.ha,craft.ca,craft.nst)
    print("Stringers (y,z) are:\n",pos)

    Zcg = calcCentroid(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst)
    print("Centroid z-coordinate is:\n", Zcg)

    Izz,Iyy = calcInertia(craft.ca,craft.ha,craft.tsk,craft.tsp,craft.tst,craft.Ast,Zcg,pos)
    print("Izz and Iyy:\n",Izz, Iyy)

    #calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,1,0)
    

    drawSection(craft.ha,craft.ca,-pos[1,:],-pos[0,:],Zcg)
if __name__ == "__main__":
    main()
