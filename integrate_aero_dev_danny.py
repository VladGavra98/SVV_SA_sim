from interpolationNew import interpolate,inner_interpolate,interpolate_data
import numpy as np
import math as m
from integration import integrationArray
import scipy.integrate
import matplotlib.pyplot as plt

zsc = -0.1115386749279281

filename = "aerodynamicloada320.dat"
data = np.loadtxt(filename, delimiter=",")
d_x = np.transpose(data)

Nx = (len(d_x[:, 0])) - 1
Nz = (len(d_x[0, :])) - 1
Ca = 0.547
la = 2.771
x = []
for i in range(Nx + 1):
    x_local = (la / 4) * ((1 - m.cos(i * m.pi / (Nx + 1))) + (1 - m.cos((i + 1) * m.pi / (Nx + 1))))
    x.append(x_local)
### Make the z-location
z = []
for i in range(Nz + 1):
    z_local = -(Ca / 4) * ((1 - m.cos(i * m.pi / (Nz + 1))) + (1 - m.cos((i + 1) * m.pi / (Nz + 1))))
    z.append(z_local)
z = z[::-1]


### Creating some new arrays to test the interpolation with
nXInterp = 50
nZInterp = 50
xgridInterp = np.linspace(x[0], x[-1], nXInterp)
zgridInterp = np.linspace(z[0], z[-1], nZInterp)
dxInterp = np.abs(x[0]-x[-1])/(nXInterp-1)
dzInterp = np.abs(z[0]-z[-1])/(nZInterp-1)
ansX,ansZ = interpolate(xgridInterp,zgridInterp)

d = ansX[3]
ax = ansX[0]
bx = ansX[1]
cx = ansX[2]
az = ansZ[0]
bz = ansZ[1]
cz = ansZ[2]


def interpolateX(x, d):
    ###initialize some things
    Nx = len(x) - 1

    a = np.zeros(len(x))
    b = np.zeros(len(x))
    c = np.zeros(len(x))

    ### Make h_i
    h = []
    for i in range(Nx):
        h.append(x[i + 1] - x[i])

    ### Initialize the lhs Solving matrix A
    A = np.zeros((Nx - 1, Nx - 1))
    for i in range(len(A[:, 0])):
        A[i, i] = ((1 / h[i]) + (1 / h[i + 1])) / 3
        if i != 0:
            A[i - 1, i] = 1 / (6 * h[i])
        if i != Nx - 2:
            A[i + 1, i] = 1 / (6 * h[i + 1])

    ### Initialize the rhs vector b
    b_rhs = np.zeros(Nx - 1)
    for i in range(Nx - 1):
        b_rhs[i] = (1 / h[i]) * (d[i + 1] - d[i]) - (1 / h[i + 1]) * (d[i + 2] - d[i + 1])

    ### Solve the system and input in Existing matrix
    M = np.insert([0.0, 0.0], 1, np.linalg.solve(A, b_rhs))
    for i in range(len(M) - 1):
        a[i] = (M[i + 1] - M[i]) / 6
        b[i] = M[i] / 2
        c[i] = (d[i + 1] - d[i]) / h[i] - (h[i] / 3) * M[i] - (h[i] / 6) * M[i + 1]

    return a,b,c


def integQZ(d):
    xVal = np.zeros(d.shape[0])
    i = 0
    for row in d:
        xVal[i] = integrationArray(row,zgridInterp[0],zgridInterp[-1],nZInterp-1)
        i += 1
    return xVal

def integQZ_Z(d):
    newD = np.zeros(d.shape)
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            newD[i,j] = d[i,j]*(zgridInterp[j]-zsc)
    xVal = np.zeros(newD.shape[0])
    k = 0
    for row in newD:
        xVal[k] = integrationArray(row, zgridInterp[0], zgridInterp[-1], nZInterp - 1)
        k += 1
    return xVal



def func(a,b,c,d,x):
    return a*x**3 + b*x**2+ c*x + d

vec = integQZ(d)
a,b,c = interpolateX(xgridInterp,vec)
print(a)
print(b)
print(c)




for i in range(len(vec)-1):
    x = np.linspace(xgridInterp[i],xgridInterp[i+1],100)
    y = a[i]*(x-xgridInterp[i])**3+b[i]*(x-xgridInterp[i])**2+c[i]*(x-xgridInterp[i])+vec[i]
    plt.plot(x,y)


plt.plot(xgridInterp,vec)
plt.show()


def integrateX(a,b,c,d,x,type):
    if type=="Single":

        for i in range(len(vec) - 1):
            x = np.linspace(xgridInterp[i], xgridInterp[i + 1], 100)
            y = a[i] * (x - xgridInterp[i]) ** 3 + b[i] * (x - xgridInterp[i]) ** 2 + c[i] * (x - xgridInterp[i]) + vec[i]
        pass
    elif type=="Single_X":
        pass
    elif type == "Double":
        pass


    elif type == "Triple_X":
        pass

    return
