# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:00:09 2020

@author: loqui
"""

import numpy as np
import math as m
from main_sim import *

#Imported geoemtry
theta = m.radians(26)
ha = 0.225 # m
tsk = 1.1/1000  # m
tsp = 2.9/1000  # m
tst = 1.2/1000  # m
hst = 15./1000   # m
wst = 20./1000   # m
nst = 17  # -
ca = 0.547
x = 2.771
x1 = 0.153  # m
x2 = 1.281  # m
x3 = 2.681  # m
xa = 0.28   # m
ha = 0.225
d1 = 0.01103  # m
d3 = 0.01642  # m
E = 72.9*10**9       # E-modulus (Pa)
G = 27.1*10**9
P = 91.7*1000

StPos = calcStPose(ha,ca,nst)
Ast = calcStArea(tst, hst, wst)
zc = calcCentroid(ha, ca, tsk, tsp, tst, hst, wst, nst)
Izz,Iyy = calcInertia(ca, ha, tsk, tsp, tst, Ast, zc, StPos)
J = 0.00024311681258111343
zsc = -0.119227056

#Shear flow? Aero load?

q = 1


#Beninf & Torsion stiffness
EIzz = E*Izz
EIyy = E*Iyy
GJ = G*J

# Line 1: Sum of forces in y-direction
line1 = np.array([1,0,np.sin(theta),1,0,1,0,0,0,0,0,0])

# Line 2: Sum of forces in the z-direction
line2 = np.array([0,1,np.cos(theta),0,1,0,1,0,0,0,0,0])

# Line 3: Moments around y
line3 = np.array([0,x-x1,np.cos(theta)*(x-(x2-xa/2)),0,x-x2,0,x-x3,0,0,0,0,0])

# Line 4: Moments around z
line4 = np.array([x-x1,0,np.sin(theta)*(x-(x2-xa/2)),x-x2,0,x-x3,0,0,0,0,0,0])

# Line 5: Torque
line5 = np.array([(-ha/2 - zsc),0,np.sin(theta)*(0 - zsc)+np.cos(theta)*ha/2,(-ha/2 - zsc),0,(-ha/2 - zsc),0,0,0,0,0,0])

# Line 6: Boundary condition 1  vy(x1) + theta(x1)*(zsc+ha/2)=d1*cos(theta0)
line6 = np.array([(1/EIzz)*(1/6)*(x1-x1)**3+(1/GJ)*(-ha/2 - zsc)*(x1-x1)*(zsc+ha/2),0,0,0,0,0,0,0,0,x1,1,(zsc+ha/2)])

# Line 7: Boundary condition 2  vy(x2) + theta(x2)*(zsc+ha/2)=0
line7 = np.array([(1/EIzz)*(1/6)*(x2-x1)**3+(1/GJ)*(-ha/2 - zsc)*(x2-x1)*(zsc+ha/2),0,(1/EIzz)*(1/6)*np.sin(theta)*(x2-(x2-xa/2))**3 + (1/GJ)*np.sin(theta)*(0 - zsc)*(x2-(x2-xa/2))*(zsc+ha/2)+(1/GJ)*np.cos(theta)*ha/2*(x2-(x2-xa/2))*(zsc+ha/2),(1/EIzz)*(1/6)*(x2-x2)**3 + (1/GJ)*(-ha/2 - zsc)*(x2-x2)*(zsc+ha/2),0,0,0,0,0,x2,1,(zsc+ha/2)])

# Line 8: Boundary condition 3  vy(x3) + theta(x3)*(zsc+ha/2)=d3*cos(theta0)
line8 = np.array([(1/EIzz)*(1/6)*(x3-x1)**3+(1/GJ)*(-ha/2 - zsc)*(x3-x1)*(zsc+ha/2),0,(1/EIzz)*(1/6)*np.sin(theta)*(x3-(x2-xa/2))**3 + (1/GJ)*np.sin(theta)*(0 - zsc)*(x3-(x2-xa/2))*(zsc+ha/2)+(1/GJ)*np.cos(theta)*ha/2*(x3-(x2-xa/2))*(zsc+ha/2),(1/EIzz)*(1/6)*(x3-x2)**3 + (1/GJ)*(-ha/2 - zsc)*(x3-x2)*(zsc+ha/2),0,(1/EIzz)*(1/6)*(x3-x3)**3 + (1/GJ)*(-ha/2 - zsc)*(x3-x3)*(zsc+ha/2),0,0,0,x3,1,(zsc+ha/2)])

# Line 9: Boundary condition 4  vz(x1) = -d1*sin(theta0)
line9 = np.array([0,(1/EIyy)*(1/6)*(x1-x1)**3,0,0,0,0,0,x1,1,0,0,0])

#Line 10: Boundary condition 5 vz(x2) = 0
line10 = np.array([0,(1/EIyy)*(1/6)*(x2-x1)**3,(1/EIyy)*np.cos(theta)*(1/6)*(x2-(x2-xa/2))**3,0,(1/EIyy)*(1/6)*(x2-x2)**3,0,0,x2,1,0,0,0])

# Line 11: Boundary condition 6 vz(x3) = -d3*sin(theta0)
line11 = np.array([0,(1/EIyy)*(1/6)*(x3-x1)**3,(1/EIyy)*np.cos(theta)*(1/6)*(x3-(x2-xa/2))**3,0,(1/EIyy)*(1/6)*(x3-x2)**3,0,(1/EIyy)*(1/6)*(x3-x3)**3,x3,1,0,0,0])

#Line 12: Boundary condition 7 for twist: Z-deformation at actuator I is 0
line12 = np.array([(1/(G*J))*(-ha/2 - zsc)*(x-x1),0,(1/(G*J))*np.sin(theta)*(0 - zsc)*(x-(x2-xa/2))+(1/(G*J))*np.cos(theta)*ha/2*(x-(x2-xa/2)),(1/(G*J))*(-ha/2 - zsc)*(x-x2),0,(1/(G*J))*(-ha/2 - zsc)*(x-x3),0,0,0,0,0,1])

#Putting together the matrix
BadBoi = np.matrix([line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12])
print(BadBoi)

B = np.matrix([[P*np.sin(theta)],
			   [P*np.cos(theta)],
			   [P*np.cos(theta)*(x-(x2+xa/2))],
			   [P*np.sin(theta)*(x-(x2+xa/2))], #-q Add torque afterwards
			   [-P*(np.cos(theta)*ha/2 + zsc*np.sin(theta))],#-q Add Torque,
			   [d1*np.cos(theta)],#+(1/(E*Izz))*q -(1/(G*J))*q Add Torque
			   [0], #[(1/(E*Izz))*q - (1/(G*J))*q] Add Torque
			   [d3*np.cos(theta) + (1/(E*Izz))*0.5*P*np.sin(theta)*(x3-(x2+xa/2))**3 -(1/(G*J))*P*(x3-(x2+xa/2))*(ha/2*np.cos(theta)+zsc*np.sin(theta)) ], # + (1/(E*Izz))*q  - (1/(G*J))*q Add Torque
			   [-d1*np.sin(theta)],
			   [0],
			   [-d3*np.sin(theta) + (1/6)*P*np.cos(theta)*(x3-(x2+xa/2))**3],
			   [(1/(G*J))*P*(x-(x2+xa/2))*(ha/2*np.cos(theta)+zsc*np.sin(theta))]]) #-(1/(G*J))*q Add Torque


Results = np.dot(np.linalg.inv(BadBoi),B)
print(Results)

#Vz
#0,(1/EIyy)*(1/6)*(x-x1)**3,(1/EIyy)*np.cos(theta)*(1/6)*(x-(x2-xa/2))**3,0,(1/EIyy)*(1/6)*(x-x2)**3,0,(1/EIyy)*(1/6)*(x-x3)**3,x,1,0,0,0;
#
#Vy
#(1/EIzz)*(1/6)*(x-x1)**3,0,(1/EIzz)*(1/6)*np.sin(theta)*(x-(x2-xa/2))**3,(1/EIzz)*(1/6)*(x-x2)**3,0,(1/EIzz)*(1/6)*(x-x3)**3,0,0,0,x,1,0;

#Theta
#(1/GJ)*(-ha/2 - zsc)*(x-x1),0,(1/GJ)*np.sin(theta)*(0 - zsc)*(x-(x2-xa/2))+(1/GJ)*np.cos(theta)*ha/2*(x-(x2-xa/2)),(1/GJ)*(-ha/2 - zsc)*(x-x2),0,(1/GJ)*(-ha/2 - zsc)*(x-x3),0,0,0,0,0,1;

#Order of the results matrix
#R=[R1y
#	R1z
#	RI
#	R2y
#	R2z
#	R3y
#	R3z
#	C1
#	C2
#	C3
#	C4
#	C5]
