# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 17:40:09 2020

@author: loqui
"""

import numpy as np

def interpolation2(n,x_vec,values):


	coeff = np.zeros([2,len(values)-1])
	h = np.diff(x_vec)


	for i in range(len(values)-1):
		coeff[0,i] = values[i]
		coeff[1,i] = (values[i+1]-values[i])/(h[i])
	return coeff


def integrate_z1(x_vec, z, abcd, x_pos):
	qx = np.array([])

# 	print(abcd)

	for i in range(x_pos+1):
		integrationz = 0
		for j in range(len(z)-1):
			integrationz += abcd[1][0][i,j]/4*(z[j+1]-z[j])**4 + abcd[1][1][i,j]/3*(z[j+1]-z[j])**3 + abcd[1][2][i,j]/2*(z[j+1]-z[j])**2 + abcd[1][3][i,j]*(z[j+1]-z[j])
		qx = np.append(qx,integrationz)
	return qx

def integrate_z2(x_vec,z,abcd,x_pos,zsc):
	qx = np.array([])
	for i in range(x_pos+1):
		integrationz = 0
		for j in range(len(z)-1):
#			integrationz += (abcd[1][0][i,j]/4*(z[j+1]-z[j])**4 + abcd[1][1][i,j]/3*(z[j+1]-z[j])**3 + abcd[1][2][i,j]/2*(z[j+1]-z[j])**2 + abcd[1][3][i,j]*(z[j+1]-z[j]))*0.5*(np.abs(zsc)+(z[j+1]-z[j]))**2
			integrationz += np.abs(zsc)*(abcd[1][0][i,j]/4*(z[j+1]-z[j])**4 + abcd[1][1][i,j]/3*(z[j+1]-z[j])**3 + abcd[1][2][i,j]/2*(z[j+1]-z[j])**2 + abcd[1][3][i,j]*(z[j+1]-z[j])) - (abcd[1][0][i,j]/5*(z[j+1]-z[j])**5 + abcd[1][1][i,j]/4*(z[j+1]-z[j])**4 + abcd[1][2][i,j]/3*(z[j+1]-z[j])**3 + abcd[1][3][i,j]/2*(z[j+1]-z[j])**2)
		qx = np.append(qx,integrationz)
	return qx



def integratefinal(x_vec,coeff,b,integration):

	h = np.diff(x_vec)
	location = 0


	for i in range(len(x_vec)):
		if x_vec[i] < b < x_vec[i+1]:
			location = i
			break
		elif x_vec[i] == b:
			location = i
			break
	if location == len(x_vec)-1:
		location = len(x_vec)-2


	x_b = b-x_vec[location]


	coeff = interpolation2(len(x_vec),x_vec,coeff)
	area_vec = np.array([0])
	totalarea = 0

	for i in range(location+1):
		if i < location:
			area = coeff[0,i]*h[i] + coeff[1,i]/2*h[i]**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)
		elif i == location:
			area = coeff[0,i]*x_b + coeff[1,i]/2*x_b**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)


	if integration==1:
		return totalarea


	coeff = interpolation2(len(x_vec),x_vec,area_vec)
	totalarea = 0
	area_vec = np.array([0])
	for i in range(location+1):
		if i < location:
			area = coeff[0,i]*h[i] + coeff[1,i]/2*h[i]**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)
		elif i == location:
			area = coeff[0,i]*x_b + coeff[1,i]/2*x_b**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)


	if integration==2:
		return totalarea


	coeff = interpolation2(len(x_vec),x_vec,area_vec)
	totalarea = 0
	area_vec = np.array([0])
	for i in range(location+1):
		if i < location:
			area = coeff[0,i]*h[i] + coeff[1,i]/2*h[i]**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)
		elif i == location:
			area = coeff[0,i]*x_b + coeff[1,i]/2*x_b**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)


	if integration==3:
		return totalarea


	coeff = interpolation2(len(x_vec),x_vec,area_vec)
	totalarea = 0
	area_vec = np.array([0])
	for i in range(location+1):
		if i < location:
			area = coeff[0,i]*h[i] + coeff[1,i]/2*h[i]**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)
		elif i == location:
			area = coeff[0,i]*x_b + coeff[1,i]/2*x_b**2
			totalarea += area
			area_vec = np.append(area_vec,totalarea)


	if integration==4:
		return totalarea

