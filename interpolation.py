# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:18:10 2020

@author: gille
"""
    


import numpy as np
import math as m
    
    
def interpolate():
    
    ### Initiating some variables and data ###
    ##################################################
    ### Import the aero data, make into 2D numpy array
    filename="aerodynamicloada320.dat"
    data = np.loadtxt(filename, delimiter=",")
    d_x = np.transpose(data)
    d_z = np.transpose(data)
    
    ### make all the interpolant data arrays
    a_x = np.zeros(np.shape(d_x))
    b_x = np.zeros(np.shape(d_x))
    c_x = np.zeros(np.shape(d_x))
    a_z = np.zeros(np.shape(d_x))
    b_z = np.zeros(np.shape(d_x))
    c_z = np.zeros(np.shape(d_x))
    
    
    ### To make code more readable...
    Nx = len(d_x[:,0])-1
    Nz = len(d_z[0,:])-1
    Ca = 0.547
    la = 2.771
    
    
    
    
    
    ### Initiating positions ####
    ##################################################
    ### Make the z-location
    z = []
    for i in range(Nz+1):
        z_local = -(Ca/4)*((1-m.cos(i*m.pi/(Nz+1)))+(1-m.cos((i+1)*m.pi/(Nz+1))))
        z.append(z_local)
        
    ### Make x-loaction
    x = []
    for i in range(Nx+1):
        x_local = (la/4)*((1-m.cos(i*m.pi/(Nx+1)))+(1-m.cos((i+1)*m.pi/(Nx+1))))
        x.append(x_local)
    
    
    
    
    ### Step 1, Get cubic splines in x-direction ###
    #################################################
    ### first make the h_i
    h = []
    for i in range(0, Nx):
        h.append(x[i+1]-x[i])
    
    for j in range(Nz+1):
        ### Initialize the lhs Solving matrix A
        A = np.zeros((Nx-1, Nx-1))
        for i in range(len(A[:,0])):
            A[i,i] = ((1/h[i])+(1/h[i+1]))/3
            if i!=0:
                A[i-1,i] = 1/(6*h[i])
            if i!=Nx-2:
                A[i+1,i] = 1/(6*h[i+1])
        
        ### Initialize the rhs vector b
        b = np.zeros(Nx-1)
        for i in range(Nx-1):
            b[i] = (1/h[i])*(d_x[i+1,j]-d_x[i,j])-(1/h[i+1])*(d_x[i+2,j]-d_x[i+1,j])
            
        ### Solve the system and input in Existing matrix
        M = np.insert([0.0,0.0], 1, np.linalg.solve(A, b))
        for i in range(len(M)-1):
            a_x[i,j] = (M[i+1]-M[i])/6
            b_x[i,j] = M[i]/2
            c_x[i,j] = (d_x[i+1,j]-d_x[i,j])/h[i] - (h[i]/3)*M[i] - (h[i]/6)*M[i+1]
    
    
    
    
    ### Step 2, get Cubic splines in z direction ###
    #################################################
    ### first make h_j
    h = []
    for j in range(0,Nz):
        h.append(z[j+1]-z[j])
        
    for i in range(Nx+1):
        ### Initialize the lhs Solving matrix A
        A = np.zeros((Nz-1, Nz-1))
        for j in range(len(A[:,0])):
            A[j,j] = ((1/h[i])+(1/h[i+1]))/3
            if j!=0:
                A[j-1,j] = 1/(6*h[i])
            if j!=Nz-2:
                A[j+1,j] = 1/(6*h[i+1])
        
        ### Initialize the rhs vector b
        b = np.zeros(Nz-1)
        for j in range(Nz-1):
            b[i] = (1/h[i])*(d_z[i,j+1]-d_z[i,j])-(1/h[i+1])*(d_z[i,j+2]-d_z[i,j+1])
            
        ### Solve the system and input in Existing matrix
        M = np.insert([0.0,0.0], 1, np.linalg.solve(A, b))
        for j in range(len(M)-1):
            a_z[i,j] = (M[j+1]-M[j])/6
            b_z[i,j] = M[j]/2
            c_z[i,j] = (d_z[i,j+1]-d_z[i,j])/h[i] - (h[i]/3)*M[i] - (h[i]/6)*M[i+1]
            
    sol_x = [a_x, b_x, c_x, d_x]
    sol_z = [a_z, b_z, c_z, d_z]
        
        
    return sol_x, sol_z


e = interpolate()

        

