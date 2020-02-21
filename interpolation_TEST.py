# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:07:36 2020

@author: gille
"""

from interpolation import interpolate
import numpy as np
import math as m
import scipy as sc
import matplotlib.pyplot as plt



### Again making some initial variable... ###
#################################################
filename="aerodynamicloada320.dat"
data = np.loadtxt(filename, delimiter=",")
d_x = np.transpose(data)  

Nx = (len(d_x[:,0]))-1
Nz = (len(d_x[0,:]))-1
Ca = 0.547
la = 2.771




### Initiating the given x- and z locations ###
##################################################
x = []
for i in range(Nx+1):
    x_local = (la/4)*((1-m.cos(i*m.pi/(Nx+1)))+(1-m.cos((i+1)*m.pi/(Nx+1))))
    x.append(x_local)
    
z = []
for i in range(Nz+1):
    z_local = (Ca/4)*((1-m.cos(i*m.pi/(Nz+1)))+(1-m.cos((i+1)*m.pi/(Nz+1))))
    z.append(z_local)

### Creating some new arrays to test the interpolation with 
x_new = np.linspace(x[0], x[-1], 100)
z_new = np.linspace(z[0], z[-1], 100)





### Doing the actual interpolation! ###
################################################
### Call the created interolation function
result = interpolate(x_new, z_new)[0][3]
### Plot it...
plt.subplot(121)
plt.title("Our cubic interpolation")
plt.imshow(result, cmap='gnuplot2')





### Check it by making scipy do the exact same kind of interpolation ###
##################################################
### First interpolate in x-direction, as we did
result_np = np.zeros((len(x_new), len(z)))
for j in range(len(z)):
    y = d_x[:,j]
    ### This specific interpolation is also 1d cubic interpolation!!
    tck = sc.interpolate.splrep(x, y, s=0)
    y_new = sc.interpolate.splev(x_new, tck, der=0)
    for i in range(len(y_new)):
        result_np[i,j] = y_new[i]

### Now interpolate in z-direction...
result_old = result_np.copy()
result_np = np.zeros((len(x_new), len(z_new)))
for i in range(len(x_new)):
    y = result_old[i,:]
    tck = sc.interpolate.splrep(z, y, s=0)
    y_new = sc.interpolate.splev(z_new, tck, der=0)
    for j in range(len(z_new)):
        result_np[i,j] = y_new[j]

### Plot results and save figure      
plt.subplot(122)
plt.title('Interpolation using Scipy')
plt.imshow(result_np, cmap='gnuplot2')
plt.savefig("cubic_vs_numpy", dpi=400)
plt.show()



### Plotting the difference
plt.imshow(result_np-result, cmap='seismic')
plt.title("Difference between interpolation methods", y=1.05)
plt.colorbar()
plt.show()




### Calculating the norm of the difference between these two calculations to check correctness
error = result.reshape((len(x_new)*len(z_new))) - result_np.reshape((len(x_new)*len(z_new)))
error = m.sqrt((np.dot(error, error)))/(len(error))
print("The average error at each node between this cubic interpolation and numpys cubic interpolation is: ", error)

