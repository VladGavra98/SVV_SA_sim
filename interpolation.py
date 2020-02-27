# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:18:10 2020
@author: gille
"""



import numpy as np
import math as m

### Inner interpolation function used often
################################################
### x will be direction of interpolation. Switch inputs if other way around ###
### Also d must then be inputted as transposed ###
def inner_interpolate(x, z, d, transpose):
    ###initialize some things
    Nx = len(x)-1
    Nz = len(z)-1

    a = np.zeros((len(x), len(z)))
    b = np.zeros((len(x), len(z)))
    c = np.zeros((len(x), len(z)))

    ### Make h_i
    h = []
    for i in range(Nx):
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
        b_rhs = np.zeros(Nx-1)
        for i in range(Nx-1):
            b_rhs[i] = (1/h[i])*(d[i+1,j]-d[i,j])-(1/h[i+1])*(d[i+2,j]-d[i+1,j])

        ### Solve the system and input in Existing matrix
        M = np.insert([0.0,0.0], 1, np.linalg.solve(A, b_rhs))
        for i in range(len(M)-1):
            a[i,j] = (M[i+1]-M[i])/6
            b[i,j] = M[i]/2
            c[i,j] = (d[i+1,j]-d[i,j])/h[i] - (h[i]/3)*M[i] - (h[i]/6)*M[i+1]

    if transpose==True:
        a = np.transpose(a)
        b = np.transpose(b)
        c = np.transpose(c)
    return a, b, c




def interpolate_data(x, x_new, abcd, j):
    return_data = []
    i2 = 0
    for i in range(len(x)-1):
        while x[i]<=x_new[i2]<=x[i+1]:
            diff = x_new[i2]-x[i]
            x_local = abcd[0][i, j]*diff**3 + abcd[1][i, j]*diff*diff + abcd[2][i, j]*diff + abcd[3][i, j]
            return_data.append(x_local)
            i2+=1
            if i2==len(x_new):
                break
    return return_data

#+++++++++++++++++++++= Main interpolation function ++++++++++++++++++++++++++++
def interpolate(x_inp, z_inp,filename):

    ### Initiating some variables and data ###
    ##################################################
    ### Import the aero data, make into 2D numpy array
    #filename="aerodynamicloada320.dat"
    data = np.loadtxt(filename, delimiter=",")
    d_x = np.transpose(data)*1000

    ### To make code more readable...
    Nx = len(d_x[:,0])-1
    Nz = len(d_x[0,:])-1
    Ca = 0.547
    la = 2.771



    ### Step 1, Initiating positions ####
    ##################################################
    ### Make x-loaction
    x = []
    for i in range(Nx+1):
        x_local = (la/4)*((1-m.cos(i*m.pi/(Nx+1)))+(1-m.cos((i+1)*m.pi/(Nx+1))))
        x.append(x_local)

    ### Make the z-location
    z = []
    for i in range(Nz+1):
        z_local = -(Ca/4)*((1-m.cos(i*m.pi/(Nz+1)))+(1-m.cos((i+1)*m.pi/(Nz+1))))
        z.append(z_local)
    z = z[::-1]

    ### A safety net to assure correct inputs ###
    if x_inp[0] < x[0] or x_inp[-1] > x[-1] or z_inp[0] < z[0] or z_inp[-1] > z[-1]:
        print("incorrect inputs given for interpolate function. Aborting...")
        return


    ### Step 2, Get cubic splines in x-direction ###
    #################################################
    a_x, b_x, c_x = inner_interpolate(x, z, d_x, False)




    ### Step 3, get data at new x points, along existing z lines ###
    ###################################################
    data_x = np.zeros((len(x_inp), Nz+1))
    for j in range(Nz+1):
        data_x[:,j] = interpolate_data(x, x_inp, [a_x, b_x, c_x, d_x], j)



    ### Step 4, create splines in z direction along new data points ###
    ##################################################
    a_z, b_z, c_z = inner_interpolate(z, x_inp, np.transpose(data_x), True)



    ### Step 5, get data points at new x AND z points
    data_z = np.zeros((len(x_inp), len(z_inp)))
    for i in range(len(x_inp)):
        data_z[i,:] = interpolate_data(z, z_inp,[np.transpose(a_z), np.transpose(b_z), np.transpose(c_z), np.transpose(data_x)], i)



    ### Step 6, create splines through these new points in both directions
    a_x, b_x, c_x = inner_interpolate(x_inp, z_inp, data_z, False)
    a_z, b_z, c_z = inner_interpolate(z_inp, x_inp, np.transpose(data_z), True)



    ### Return answer ###
    ans_x = [a_x, b_x, c_x, data_z]
    ans_z = [a_z, b_z, c_z, data_z]
    return ans_x, ans_z