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

def calcCircum(ha,ca):
    #Stick to the name given in the flow chart OR
    # cleary write what you cahnged
    """Input: ha, ca
       Output: lCirc = length of the section circumference
    """
    return np.pi*ha/2 + 2*(np.sqrt((ca-ha/2)**2 + (ha/2)**2))

class Aircraft:
    def __init__(self,name):
        if name=="A320" or name=="a320":
            self.la = 2.71 #m
            self.ca = 0.5
            self.ha = 0.2
            self.tsk = 0.00011
            self.tsp = 0.0002
            self.wst = 0.002
            self.hst = 0.002
            self.nst = 17

def integration(function,n,a,b):
	
	z0 = a 											# Lower boundary of integration
	zf = b 											# Upper boundary of integration
	n = n 											# Number of steps
	deltaz = round((zf-z0)/n,8) 					# Length of each step
	zvector = np.array([0])							# Summation of intervals
	h = deltaz/2 									# Half the step
	primvector = np.array([0]) 						# Results from Simpson's rule
	valuez = 0
	
	for i in range(n): 								# Calcualtion of summation of intervals
		valuez += deltaz
		zvector = np.append(zvector,valuez)
	
			
	for j in range(n):								# Simpson's rule
		value = h/3*(function(zvector[j]) + 4*function(zvector[j]+h) + function(zvector[j+1]))
		primvector = np.append(primvector,value)
	
	A = np.empty([n+1,n+1])							# Interpolation
	for i in range(n+1):
		for j in range(n+1):
			A[i,j] = zvector[i]**j
		
	weights = np.dot(np.linalg.inv(A),np.transpose(primvector)) 	# Weights of interpolation, value in front of the xs
	total = np.sum(primvector)										# Total value of integration
	
	
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
#++++++++++++++++++++++++++++ Main +++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
    aircraft = Aircraft("A320")
    print("Circumference: \n",calcCircum(aircraft.ha,aircraft.ca))
    drawSection(aircraft.ha,aircraft.ca)


if __name__ == "__main__":
    main()
