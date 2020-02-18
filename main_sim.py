# -*- coding: utf-8 -*-
"""
SVV 2020- Structural Analysis Assignment

Simulation for stress and deflection in A320 Airleron

@author: vladg
"""
import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
plt.close('all')

#++++++++++++++++++++++++++ Constants +++++++++++++++++++++++++++++++++++++++++++++++++++
g = 9.81 #m/s2


def calcShFlow(ha,ca,tsk,tsp, tst, hst, wst,nst,Sz,Sy):
    #open section shear flow

    # testnumbers
    Izz = 1
    Iyy = 1
    ca = 7
    ha = 6

    zCentroid = calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst)
    stArea = calcStArea(tst,hst,wst)
    stringerPos = calcStPose(ha, ca, nst)
    ##### top flat plate
    # filter stringers on top flat plate


    plateYLength = ha/2
    plateZLength = ca-ha/2
    plateLength = np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))

    ############ (1) top flate plate starts here
    n1 = 100
    stringerPosFilt1 = stringerPos[:, (stringerPos[0, :] >= 0) & (stringerPos[1, :] <= -ha / 2)]
    sVec1 = np.linspace(0,plateLength,n1+1)
    yVec1 = (plateYLength/plateLength) * sVec1
    zVec1 = -zCentroid -ca + (plateZLength/plateLength) * sVec1
    ds1 = plateLength/n1

    intYVec1,intZVec1 = calcIntegralArray(zVec1,yVec1,sVec1,n1)

    #intYVec1 = addStringerContribution(intYVec1,yVec1,zVec1,ds1,n1,stringerPosFilt1,stArea,"Y")
    #intZVec1 = addStringerContribution(intZVec1,yVec1,zVec1,ds1,n1,stringerPosFilt1,stArea,"Z")
    qs1 = -(Sz/Iyy)*intZVec1 - (Sy/Izz)*intYVec1
    print(qs1)

    ##########  (2) spar plate starts here
    n2 = 100
    #no stringers on spar duh
    sVec2 = np.linspace(0,ha,n1+1)
    yVec2 = ha/2-sVec2
    zVec2 = np.ones(len(sVec2))*(-zCentroid -ha/2)
    ds2 = ha/n2

    intYVec2,intZVec2 = calcIntegralArray(zVec2,yVec2,sVec2,n2)
    qs2 = -(Sz/Iyy)*intZVec2 - (Sy/Izz)*intYVec2
    #add last value from qs1
    qs2 += qs1[-1]
    print(qs2)

    ######## (3) lower flat plate starts here
    n3 = 100
    stringerPosFilt3 = stringerPos[:, (stringerPos[0, :] <= 0) & (stringerPos[1, :] <= -ha / 2)]
    sVec3 = np.linspace(0, plateLength, n3 + 1)
    yVec3 = (-ha/2) + (plateYLength / plateLength) * sVec3
    zVec3 = (-zCentroid-ha/2) - (plateZLength / plateLength) * sVec3
    ds3 = plateLength / n3

    intYVec3, intZVec3 = calcIntegralArray(zVec3, yVec3, sVec3, n3)

    #intYVec3 = addStringerContribution(intYVec3, yVec3, zVec3, ds3, n3, stringerPosFilt3, stArea, "Y")
    #intZVec3 = addStringerContribution(intZVec3, yVec3, zVec3, ds3, n3, stringerPosFilt3, stArea, "Z")
    qs3 = -(Sz / Iyy) * intZVec3 - (Sy / Izz) * intYVec3
    #add last value from qs2
    qs3 += qs2[-1]

    print(qs3)

    return



def calcIntegralArray(z,y,s, n):
    intYVec = np.array([0])
    intZVec = np.array([0])
    for i in range(1, n):
        integralY = sp.integrate.simps(y[:i + 1], s[:i + 1])
        intYVec = np.append(intYVec, integralY)
        integralZ = sp.integrate.simps(z[:i + 1], s[:i + 1])
        intZVec = np.append(intZVec, integralZ)
    return intYVec,intZVec


def addStringerContribution(integrated,yVec,zVec,ds,n,stringerPos,stArea,direction):
    used = np.array([])
    newIntegrated = integrated
    for i in range(0,n+1):
        for j in range(0,stringerPos.shape[1]):
            dist = calcDist(stringerPos[0,j],stringerPos[1,j],yVec[i],zVec[i])
            if dist < ds and j not in used:
                used = np.append(used,j)
                if direction == "Y":
                    newIntegrated[j+1:] += stArea*stringerPos[0,j]
                elif direction == "Z":
                    newIntegrated[j + 1:] += stArea * stringerPos[1, j]
    if len(used) != stringerPos.shape[1]:
        print("Warning addStringerContribution(): Not all stringers have been used")
        print(used)
    return newIntegrated

def calcDist(y1,z1,y2,z2):
    return np.sqrt(np.power(y1-y2, 2) + np.power(z1-z2, 2))


def drawSection(ha,ca):
    """Plots the cross-section."""
    fig,ax = plt.subplots()
    plt.title("Cross section")
    lines = [[(0, 0), (ha/2, ha/2)], [(0, 0), (ha/2, -ha/2)], [(ha/2, ha/2), (ha/2, -ha/2)],
             [(ha/2,ha/2),(ca,0)],[(ha/2,-ha/2),(ca,0)]]

    lc = mc.LineCollection(lines, linewidths=2)
    circle = plt.Circle((ha/2,0),ha/2,color='b',fill=False)
    ax.add_artist(circle)
    ax.add_collection(lc)
    ax.autoscale()

def calcStPose(ha,ca,nst):
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
        pos[:,int(nst/2)+1 - i - 1] = [ha/2 * np.sin(np.radians(i*unit/(ha/2))),-ha/2 * np.cos(np.radians(i*unit/(ha/2)))]
        pos[:,int(nst/2)+1 + i -1 ] = [-ha/2 * np.sin(np.radians(i*unit/(ha/2))),-ha/2 * np.cos(np.radians(i*unit/(ha/2)))]
     #Going along the skin
    if (2*stCircle +1 < nst):
        left   = unit - (ha/2 * np.pi/2 - stCircle*unit)    #left distance from the unit
        leftSt = int(nst/2) - stCircle     #nr of stringers left to palce outside the circle
        alfa   = np.arctan2(ha/2,ca-ha/2)   #the slope angle of the straight part
        lineSt = int((semiCircum - (ha/2 * np.pi/2) - left)/leftSt)
        leftTE = semiCircum- (ha/2 * np.pi/2)- lineSt*unit
        for i in range(leftSt):
             pos[:,i]  = [(leftTE+unit*i)*np.sin(alfa),-ca + (leftTE+unit*i)*np.sin(alfa)]
             pos[:,len(pos[1,:])-i-1] = [-(leftTE+unit*i)*np.sin(alfa),-ca + (leftTE+unit*i)*np.sin(alfa)]
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


def calcCentroid(ha,ca,tsk,tsp,tst,hst,wst,nst):
    stArea = calcStArea(tst,hst,wst)
    stPos = calcStPose(ha,ca,nst)
    plateYLength = ha / 2
    plateZLength = ca - ha / 2
    sumStAreaZ = np.sum(stPos[1,:]*stArea)
    sumAreaZ = np.pi*tsk*(ha/2) * (-ha/2+(2/np.pi)*(ha/2)) + np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))*tsk*2 * (-ha/2 - plateZLength/2) + ha*tsp*(-ha/2) + sumStAreaZ
    sumArea = np.pi*tsk*(ha/2) + np.sqrt(np.power(plateYLength, 2) + np.power(plateZLength, 2))*tsk*2 + ha*tsp + stArea*nst
    zCentroid = sumAreaZ/sumArea
    return zCentroid

def calcInertia(Ca,H,Tsk,Tsp,Tst,Ast,Zcg,StPos):
    
    #Length of the straight skin
    Lsk = np.linalg.norm([H/2,Ca-H/2])
    
    
    # ------------------------   Izz   ------------------------
    # I_zz consists of 4 parts: spar (1), skin plates (2), skin semicircular (3), stiffeners (4)
    
    I_zz1 = 1/12*Tsp*H**3
    Beta_plate = atan((H/2)/(Ca-H/2))
    I_zz2 = 1/12*Tsk*(2*Lsk)**3*(sin(Beta_plate))**2
    I_zz3 = 1/128*pi*Tsk**4 
    I_zz4 = Ast*sum(StPos[0,:]**2) #calcStPos gives list of coordinates (y,z)
    
    Izz = I_zz1+I_zz2+I_zz3+I_zz4
    
    # ------------------------   Iyy   -------------------------
    # I_yy consists of 4 parts: skin plates (1), skin semicircular (2), stiffeners (3), spar (4)
    # the MoI of the thinwalled semicircle about diameter is r^3*t*pi/4 (calculated by hand)
    # the MoI of the thinwalled semicircle about cg is r^3*t*(pi/4 - 4/pi) (calculated by hand)
    
    I_yy_plate = 1/12*Tsk*(Lsk)**3*(cos(Beta_plate))**2 + Lsk*Tsk*(H/2+0.5*(Ca-H/2)-Zcg)**2
    I_yy1 = 2*I_yy_plate
    I_yy2 = (H/2)**3*Tsk*(pi/4 - 4/pi) + pi*H/2*Tsk * ((H/2-H/pi)-Zcg)**2
    I_yy3 = Ast*sum((z-Zcg)**2 for z in StPos[1,:])   #calcStPos gives list of coordinates (y,z)
    I_yy4 = H*Tsp*(H/2-Zcg)**2                      #only steiner term due to thin walled approx
    
    Iyy = I_yy1+I_yy2+I_yy3+I_yy4
        
    return Izz, Iyy

  
  #+++++++++++++++++++++++++++++++++++++ Numerical Integration ++++++++++++++++++++++++++++++++++++++++++++++++++
def integration(function,n,a,b):
	z0 = a
	zf = b
	n = n
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



def calcStArea(Tst, Hst, Wst): 
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
            self.wst = 0.02
            self.hst = 0.015
            self.nst = 17
            self.tsp = 0.0029
            self.theta = np.radians(26)  #rad
            
#++++++++++++++++++++++++++++ Main +++++++++++++++++++++++++++++++++++++++++++++++++++

def main():
    np.set_printoptions(precision=3)
    craft = Aircraft("A320")
    print("Circumference: \n",calcCircum(craft.ha,craft.ca))
    print("Stringer positions are:\n",calcStPose(craft.ha,craft.ca,17))
    print("Centroid z-coordinate is:\n", calcCentroid(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst))
    calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,1,0)

if __name__ == "__main__":
    main()
    

