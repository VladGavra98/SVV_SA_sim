# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:32:04 2020
Integration function
@author: Luis, dannywhuang
@version: 21-02
"""

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def f(x):
    return x*x


def integration(function,a,b,n):
    """Input: f,a,b, n
     Output: wights"""
    z0 = a 											# Lower boundary of integration
    zf = b 											# Upper boundary of integration

    if z0>zf:
        aux = zf
        zf  = z0
        z0  = aux

    deltaz     = round((zf-z0)/n,8) 					# Length of each step
    zvector    = np.linspace(z0,zf,n+1)				# Summation of interval
    primvector = np.zeros((len(zvector))) 		     # Results from Simpson's rule
    # h = deltaz/2       # Half the step
    valuez = 0
    h = deltaz/2

    #  	for i in range(n): 								# Calcualtion of summation of intervals
    # 		valuez += deltaz
    # 		zvector = np.append(zvector,valuez)


    for j in range(1,n):								# Simpson's rule
    	value = h/3*(function(zvector[j]) + 4*function(zvector[j]+h) + function(zvector[j+1]))
    	primvector[j] = primvector[j-1] + value

    A = np.empty([n+1,n+1])							# Interpolation
    for i in range(n+1):
    	for j in range(n+1):
    		A[i,j] = zvector[i]**j


    weights = np.dot(np.linalg.inv(A),np.transpose(primvector)) 	# Weights of interpolation, value in front of the xs
    total = primvector[-2]										# Total value of integration


    	####  Verification ######
    	#vervec = np.array([0])
    	#verificationtotal = scp.integrate.simps(function(zvector),zvector)
    	#for i in range(n):
    		#vector = np.array([zvector[i],zvector[i+1]])
    		#vervec = np.append(vervec,scp.integrate.simps(function(vector),vector))

    	#plt.plot(zvector,vervec)
    	#plt.plot(zvector,primvector)
    	#plt.show()

    return primvector[:-1],total

def integrationArray(y,a,b,n):
    """Integration with arrays as input"""
    # y = array with function values
    # a = lower bound of integration
    # b = upper bound of integration
    # n = amount of grids, n + 1 = amount of grid points

    Ny = len(y)

    if Ny != n+1:
        raise ValueError("Array length is not the same as the amount of grid points")

    if b>=a:
        h = (b-a)/n
        z0 = a
        zf = b
    elif b<a:
        h = (a-b)/n
        z0 = b
        zf = a

    if Ny%2==0:
        # if number of points is even, take average of
        # trapezoidal on first two points, simpson's on the rest
        # trapezoidal on last two points, simpson's on the rest
        value1 = (h/2)*(y[1]+y[0])
        value1 += (h/3)*np.sum(y[1:-2:2]+4*y[2:-1:2]+y[3::2])

        value2 = (h/2)*(y[-1]+y[-2])
        value2 += (h/3)*np.sum(y[0:-3:2]+4*y[1:-2:2]+y[2:-1:2])

        return (value1+value2)/2

    elif Ny%2==1:
        # if number is uneven just simpson's rule is used

        value = (h/3)*np.sum(y[0:-1:2]+4*y[1:-1:2]+y[2::2])
        return value

#Verification of integrationArra
#x = np.linspace(0,1,2)
#y = x*x
#
#print(integrationArray(y,0,1,1))
#print(scipy.integrate.simps(y,x))