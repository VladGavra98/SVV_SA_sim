# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:32:04 2020

@author: Luis
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as scp

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
  
