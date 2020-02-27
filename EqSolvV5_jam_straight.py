import numpy as np
import math as m
import matplotlib.pyplot as plt
import Final_Integration as FI
from interpolationValidation import interpolate
from validation import deflectionAlongX

Ca = 0.605
la = 2.661
x = 2.659
x1 = 0.172  # m
x2 = 1.211  # m
x3 = 2.591  # m
xa = 0.35  # m
ha = 0.205
d1 = 0  # m
d3 = 0  # m
theta = m.radians(28)
E = 73.1 * 10 ** 9  # E-modulus (Pa)
G = 28 * 10 ** 9
P = 97.4 * 1000
xI = x2 - xa / 2
xII = x2 + xa / 2
Izz = 1.0280189203385743e-05
Iyy = 8.651211860639687e-05
J = 1.5101498390705808e-05
zsc = -0.10191789189921291

# x = 2.691  # m
# x1 = 0.174  # m
# x2 = 1.051  # m
# x3 = 2.512  # m
# xa = 0.30   # m
# ha = 0.248  # m
# tsk = 1.1/1000  # m
# tsp = 2.2/1000  # m
# tst = 1.2/1000  # m
# hst = 15/100   # m
# wst = 30/100  # m
# nst = 11  # -
# d1 = 1.034/100  # m
# d3 = 2.066/100  # m
# theta = m.radians(25)  # rad
# P = 20.6*1000  #
# E = 72.9*10**9       # E-modulus (Pa)
# G = 27.1*10**9       # G-modulus (Pa)
# Iyy = 0.0001888816058419491
# Izz = 4.577138260372715e-05
# J = 1.9193311985303668e-05
# zsc = -0.12532244279667276
# xI = x2-xa/2
# xII = x2 + xa/2

# Ca = 0.505 #m
# x = 1.611 #m
# x1 = 0.125 #m
# x2 = 0.498 #m
# x3 = 1.494 #m
# xa = 0.245 #m
# ha = 0.161 #m
# d1 = 0.00389 # m
# d3 = 0.01245  # m
# P = 49.2*1000  # N
# theta = m.radians(30) #rad
# zsc = -0.08553893540215983 #m
# Izz = 4.753851442684436e-06
# Iyy = 4.5943507864451845e-05
# J = 7.748548555816593e-06
# E = 72.9*10**9
# G = 27.1*10**9
# xI = x2 - xa/2
# xII = x2 + xa/2


# Order of the results matrix
# R=[R1y 0
#	R1z 1
#	RI  2
#	R2y 3
#	R2z 4
#	R3y 5
#	R3z 6
#	C1  7
#	C2  8
#	C3  9
#	C4  10
#	C5] 11

### Again making some initial variable... ###
#################################################
filename = "aerodynamicloada320.dat"
data = np.loadtxt(filename, delimiter=",")
d_x = np.transpose(data)

Nx = (len(d_x[:, 0])) - 1
Nz = (len(d_x[0, :])) - 1

### Initiating the given x- and z locations ###
##################################################
x_vec = []
for i in range(Nx + 1):
    x_local = (la / 4) * ((1 - m.cos(i * m.pi / (Nx + 1))) + (1 - m.cos((i + 1) * m.pi / (Nx + 1))))
    x_vec.append(x_local)

z_vec = []
for i in range(Nz + 1):
    z_local = -(Ca / 4) * ((1 - m.cos(i * m.pi / (Nz + 1))) + (1 - m.cos((i + 1) * m.pi / (Nz + 1))))
    z_vec.append(z_local)
z_vec = z_vec[::-1]

## Creating some new arrays to test the interpolation with

abcd = interpolate(x_vec, z_vec)
qx = FI.integrate_z1(x_vec, z_vec, abcd, len(x_vec) - 1)
qx2 = FI.integrate_z2(x_vec, z_vec, abcd, len(x_vec) - 1, zsc)
qsy = FI.integratefinal(x_vec, qx, x_vec[-1], 1)
qmz = FI.integratefinal(x_vec, qx, x_vec[-1], 2)
qt = FI.integratefinal(x_vec, qx2, x_vec[-1], 1)
qvy1 = FI.integratefinal(x_vec, qx, x1, 4)
qvy2 = FI.integratefinal(x_vec, qx, x2, 4)
qvyII = FI.integratefinal(x_vec, qx, xII, 4)
qvy3 = FI.integratefinal(x_vec, qx, x3, 4)
qtheta1 = FI.integratefinal(x_vec, qx2, x1, 2)
qtheta2 = FI.integratefinal(x_vec, qx2, x2, 2)
qthetaII = FI.integratefinal(x_vec, qx2, xII, 2)
qtheta3 = FI.integratefinal(x_vec, qx2, x3, 2)

## Line 1: Sum of forces in y-direction
line1 = np.array([1, 0, -np.sin(theta), 1, 0, 1, 0, 0, 0, 0, 0, 0])

## Line 2: Sum of forces in the z-direction
line2 = np.array([0, 1, -np.cos(theta), 0, 1, 0, 1, 0, 0, 0, 0, 0])

## Line 3: Moments around y
line3 = np.array([0, (x - x1), -np.cos(theta) * (x - (xI)), 0, (x - x2), 0, (x - x3), 0, 0, 0, 0, 0])

## Line 4: Moments around z
line4 = np.array([x - x1, 0, -np.sin(theta) * (x - (xI)), x - x2, 0, x - x3, 0, 0, 0, 0, 0, 0])

## Line 5: Torque
line5 = np.array(
    [-(np.abs(zsc) - ha / 2), 0, np.sin(theta) * (np.abs(zsc)) - np.cos(theta) * ha / 2, -(np.abs(zsc) - ha / 2), 0,
     -(np.abs(zsc) - ha / 2), 0, 0, 0, 0, 0, 0])

## Line 6: Boundary condition 1  vy(x1) + theta(x1)*(zsc+ha/2)=d1*cos(theta0)
line6 = np.array(
    [-(1 / (6 * E * Izz)) * (x1 - x1) ** 3 - (1 / (G * J)) * (zsc + ha / 2) * (x1 - x1) * (np.abs(zsc) - ha / 2), 0, 0,
     0, 0, 0, 0, 0, 0, x1, 1, (zsc + ha / 2)])

## Line 7: Boundary condition 2  vy(x2) + theta(x2)*(zsc+ha/2)=0
line7 = np.array(
    [-(1 / (6 * E * Izz)) * (x2 - x1) ** 3 - (1 / (G * J)) * (zsc + ha / 2) * (x2 - x1) * (np.abs(zsc) - ha / 2), 0,
     -(1 / (6 * E * Izz)) * -np.sin(theta) * (x2 - (xI)) ** 3 + (1 / (G * J)) * (
                 np.sin(theta) * np.abs(zsc) * (x2 - (xI)) - np.cos(theta) * ha / 2 * (x2 - (xI))) * (zsc + ha / 2),
     -(1 / (6 * E * Izz)) * (x2 - x2) ** 3 - (1 / (G * J)) * (np.abs(zsc) - ha / 2) * (x2 - x2) * (zsc + ha / 2), 0, 0,
     0, 0, 0, x2, 1, (zsc + ha / 2)])

## Line 8: Boundary condition 3  vy(x3) + theta(x3)*(zsc+ha/2)=d3*cos(theta0)
line8 = np.array(
    [-(1 / (6 * E * Izz)) * (x3 - x1) ** 3 - (1 / (G * J)) * (zsc + ha / 2) * (x3 - x1) * (np.abs(zsc) - ha / 2), 0,
     -(1 / (6 * E * Izz)) * -np.sin(theta) * (x3 - (xI)) ** 3 + (1 / (G * J)) * (
                 np.sin(theta) * np.abs(zsc) * (x3 - (xI)) - np.cos(theta) * ha / 2 * (x3 - (xI))) * (zsc + ha / 2),
     -(1 / (6 * E * Izz)) * (x3 - x2) ** 3 - (1 / (G * J)) * (np.abs(zsc) - ha / 2) * (x3 - x2) * (zsc + ha / 2), 0,
     -(1 / (6 * E * Izz)) * (x3 - x3) ** 3 - (1 / (G * J)) * (zsc + ha / 2) * (x3 - x3) * (np.abs(zsc) - ha / 2), 0, 0,
     0, x3, 1, (zsc + ha / 2)])

## Line 9: Boundary condition 4  vz(x1) = -d1*sin(theta0)
line9 = np.array([0, -(1 / (6 * E * Iyy)) * (x1 - x1) ** 3, 0, 0, 0, 0, 0, x1, 1, 0, 0, 0])

## Line 10: Boundary condition 5 vz(x2) = 0
line10 = np.array([0, -(1 / (6 * E * Iyy)) * (x2 - x1) ** 3, (1 / (6 * E * Iyy)) * np.cos(theta) * (x2 - (xI)) ** 3, 0,
                   -(1 / (6 * E * Iyy)) * (x2 - x2) ** 3, 0, 0, x2, 1, 0, 0, 0])

## Line 11: Boundary condition 6 vz(x3) = -d3*sin(theta0)
line11 = np.array([0, -(1 / (6 * E * Iyy)) * (x3 - x1) ** 3, (1 / (6 * E * Iyy)) * np.cos(theta) * (x3 - (xI)) ** 3, 0,
                   -(1 / (6 * E * Iyy)) * (x3 - x2) ** 3, 0, -(1 / (6 * E * Iyy)) * (x3 - x3) ** 3, x3, 1, 0, 0, 0])

## Line 12: Boundary condition 7 vz(xI)*cos(theta) + vy(xI)*sin(theta) + theta(xI)*abs(zsc)*sin(theta) = 0
line12 = np.array([-(1 / (6 * E * Izz)) * np.sin(theta) * (xI - x1) ** 3 - (1 / (G * J)) * np.sin(theta) * np.abs(
    zsc) * (np.abs(zsc) - ha / 2), -(1 / (6 * E * Iyy)) * np.cos(theta) * (xI - x1) ** 3,
                   -(1 / (6 * E * Iyy)) * np.cos(theta) * np.cos(theta) * (xI - xI) ** 3 - (1 / (6 * E * Izz)) * np.sin(
                       theta) * np.sin(theta) * (xI - xI) ** 3 + (1 / (G * J)) * np.abs(zsc) * np.sin(theta) * np.sin(
                       theta) * np.abs(zsc) * (xI - xI), 0, 0, 0, 0, xI * np.cos(theta), np.cos(theta),
                   xI * np.sin(theta), np.sin(theta), -np.abs(zsc) * np.sin(theta)])

BadBoi = np.matrix([line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12])

B = np.matrix([[-P * np.sin(theta) - qsy],
			   [-P * np.cos(theta)],
			   [-P * np.cos(theta) * (x - (xII))],
			   [-P * np.sin(theta) * (x - (xII)) - qmz],  # -q Add Aero afterwards
			   [-P * (np.cos(theta) * ha / 2 - np.abs(zsc) * np.sin(theta)) - qt],  # -q Add Aero,
			   [d1 * np.cos(theta) + (1 / (E * Izz)) * qvy1 - (1 / (G * J)) * qtheta1 * (zsc + ha / 2)],
			   # +(1/(E*Izz))*q -(1/(G*J))*q Add Aero
			   [(1 / (E * Izz)) * qvy2 - (1 / (G * J)) * qtheta2 * (zsc + ha / 2)],
			   # [(1/(E*Izz))*q - (1/(G*J))*q] Add Aero
			   [d3 * np.cos(theta) + (1 / (6 * E * Izz)) * P * np.sin(theta) * (x3 - (xII)) ** 3 - (1 / (G * J)) * P * (
                           zsc + ha / 2) * (x3 - (xII)) * (np.abs(zsc) * np.sin(theta) - ha / 2 * np.cos(theta)) + (
                            1 / (E * Izz)) * qvy3 - (1 / (G * J)) * qtheta3 * (zsc + ha / 2)],
			   # + (1/(E*Izz))*q  - (1/(G*J))*q Add Aero
			   [-d1 * np.sin(theta)],
			   [0],
			   [-d3 * np.sin(theta) + (1 / (6 * E * Iyy)) * P * np.cos(theta) * (x3 - (xII)) ** 3],
			   [(1 / (E * Izz)) * qvyII * np.sin(theta) + (1 / (G * J)) * qthetaII * np.sin(theta) * np.abs(
                   zsc)]])  # -(1/(G*J))*q Add Aero

R = np.dot(np.linalg.inv(BadBoi), B)
R1 = np.linalg.solve(BadBoi, B)


def positivo(x, power):
    if power > 0 and x >= 0:
        return x ** power
    elif power == 0 and x >= 0:
        return 1
    else:
        return 0


def Sy(x):
    return R[0] * positivo(x - x1, 0) - R[2] * np.sin(theta) * positivo(x - (xI), 0) + R[3] * positivo(x - x2,
                                                                                                       0) + P * np.sin(
        theta) * positivo(x - (xII), 0) + R[5] * positivo(x - x3, 0) + FI.integratefinal(x_vec, qx, x, 1)


def Sz(x):
    return R[1] * positivo(x - x1, 0) - R[2] * np.cos(theta) * positivo(x - (xI), 0) + R[4] * positivo(x - x2,
                                                                                                       0) + P * np.cos(
        theta) * positivo(x - (xII), 0) + R[6] * positivo(x - x3, 0)


def My(x):
    return R[1] * positivo(x - x1, 1) - R[2] * np.cos(theta) * positivo(x - (xI), 1) + R[4] * positivo(x - x2,
                                                                                                       1) + P * np.cos(
        theta) * positivo(x - (xII), 1) + R[6] * positivo(x - x3, 1)


def Mz(x):
    return R[0] * positivo(x - x1, 1) - R[2] * np.sin(theta) * positivo(x - (xI), 1) + R[3] * positivo(x - x2,
                                                                                                       1) + P * np.sin(
        theta) * positivo(x - (xII), 1) + R[5] * positivo(x - x3, 1) + FI.integratefinal(x_vec, qx, x, 2)


def T(x):
    return -R[0] * (np.abs(zsc) - ha / 2) * positivo(x - x1, 0) + R[2] * (
                np.sin(theta) * np.abs(zsc) - np.cos(theta) * ha / 2) * positivo(x - (xI), 0) - R[3] * (
                       np.abs(zsc) - ha / 2) * positivo(x - x2, 0) - P * (
                       np.sin(theta) * np.abs(zsc) - np.cos(theta) * ha / 2) * positivo(x - (xII), 0) - R[5] * (
                       np.abs(zsc) - ha / 2) * positivo(x - x3, 0) + FI.integratefinal(x_vec, qx2, x, 1)


def Vy(x):
    return -(1 / (6 * E * Izz)) * (
                R[0] * positivo(x - x1, 3) - R[2] * np.sin(theta) * positivo(x - (xI), 3) + R[3] * positivo(x - x2,
                                                                                                            3) + P * np.sin(
            theta) * positivo(x - (xII), 3) + R[5] * positivo(x - x3, 3)) + R[9] * x + R[10] - (
                       1 / (E * Izz)) * FI.integratefinal(x_vec, qx, x, 4)


def Vz(x):
    return -(1 / (6 * E * Iyy)) * (
                R[1] * positivo(x - x1, 3) - R[2] * np.cos(theta) * positivo(x - (xI), 3) + R[4] * positivo(x - x2,
                                                                                                            3) + P * np.cos(
            theta) * positivo(x - (xII), 3) + R[6] * positivo(x - x3, 3)) + R[7] * x + R[8]


def Theta(x):
    return (1 / (G * J)) * (-R[0] * (np.abs(zsc) - ha / 2) * positivo(x - x1, 1) + R[2] * (
                np.sin(theta) * np.abs(zsc) - np.cos(theta) * ha / 2) * positivo(x - (xI), 1) - R[3] * (
                                        np.abs(zsc) - ha / 2) * positivo(x - x2, 1) - P * (
                                        np.sin(theta) * np.abs(zsc) - np.cos(theta) * ha / 2) * positivo(x - (xII), 1) -
                            R[5] * (np.abs(zsc) - ha / 2) * positivo(x - x3, 1)) - R[11] + (
                       1 / (G * J)) * FI.integratefinal(x_vec, qx2, x, 2)


xcheck = np.linspace(0, x, 100)

DeflectionY = np.array([])
DeflectionZ = np.array([])
for i in range(len(xcheck)):
    DeflectionY = np.append(DeflectionY, Vy(xcheck[i]) + Theta(xcheck[i]) * (ha / 2 + zsc))
    DeflectionZ = np.append(DeflectionZ, Vz(xcheck[i]))

## Check list
# Torque is good, the rest not so much, why???
# Check the Twist equation, take constant out and multiply by 2 and you get a close result\
deflectionAlongX(0, 0, xcheck, DeflectionY, DeflectionZ, "jam_straight")

## Check list
# Torque is good, the rest not so much, why???
# Check the Twist equation, take constant out and multiply by 2 and you get a close result