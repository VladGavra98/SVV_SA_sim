# -*- coding: utf-8 -*-
"""
SVV 2020- Structural Analysis Assignment

Simulation for stress and deflection in A320 Airleron


@author: vladg
@version: 21-02-#2

"""
import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt
from shear_flow import *
from matplotlib import collections  as mc
plt.close('all')
np.set_printoptions(precision=3)

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
        lineSt = int((semiCircum - (ha/2 * np.pi/2) - left)/unit)

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

def VonMisses(sigma,tau):
    """Returns the combined loading stress according to von Misses yeild criterion.
    Input:
           sigma = [[z],[y],[sigma_x]]

            """
    A = np.power(sigma[2,:],2)
    B = np.power(tau,2)
    return np.sqrt(A/2 + 3*B)

def calcStArea(Tst, Hst, Wst):   #Verified by Vlad!
    #Calculates area of stringer in m^2
    StArea = Tst * (Hst + Wst)
    return StArea

#++++++++++++++++++++++++++++++ Aircraft Class ++++++++++++++++++++++++++++++++++++++++++++
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

class Discretization:
    def __init__(self):
        self.n1 = 1000
        self.n2 = 1000
        self.n3 = 1000
        self.n4 = 1000

#++++++++++++++++++++++++++++ Main +++++++++++++++++++++++++++++++++++++++++++++++++++

def main():

    craft = Aircraft("A320")
    discret = Discretization()

    #print("Circumference: \n",calcCircum(craft.ha,craft.ca))

    stArea = calcStArea(craft.tst,craft.hst,craft.wst)
    #print("Stringer Area is:\n",stArea)

    pos = calcStPose(craft.ha,craft.ca,craft.nst)
    #print("Stringers (y,z) are:\n",pos)

    Zcg = calcCentroid(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst)
    #print("Centroid z-coordinate is:\n", Zcg)

    Izz,Iyy = calcInertia(craft.ca,craft.ha,craft.tsk,craft.tsp,craft.tst,craft.Ast,Zcg,pos)
    #print("Izz and Iyy:\n",Izz, Iyy)

    q = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,1,0, discret.n1,discret.n2,discret.n3,discret.n4)
    print("Shear flows are:\n", q)

    Zsc = calcShCenter(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,discret.n1,discret.n2,discret.n3,discret.n4)
    #print("Shear center z-coordinate is:\n", zShear)

    vm  = VonMisses(np.array([[0],[0],[0]]),q)
    print("Von Misses stress are:\n",vm)

    drawSection(craft.ha,craft.ca,-pos[1,:],-pos[0,:],Zcg,Zsc)

if __name__ == "__main__":
    main()
    plt.show()
