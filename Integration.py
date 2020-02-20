# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:32:04 2020
Integgration function
@author: Luis
@version: 1 verfied
"""
import numpy as np
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
