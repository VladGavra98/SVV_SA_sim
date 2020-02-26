import numpy as np
from main_sim import *
import matplotlib.pyplot as plt
import pandas as pd

def shearFlowGraph1(craft):
    # draw graph for Sy = 1, Sz = 0 to compare with verification model
    #verification model shear flow, Sy = 1, Sz = 0
    data = pd.read_csv("verification_data/SyOne.csv",index_col=0)
    q1Ver = data.loc[:,'q1f']
    q2Ver= data.loc[:,'q2f']
    q3Ver= data.loc[:,'q3f']
    q4Ver= data.loc[:,'q4f']
    q5Ver= data.loc[:,'q5f']
    q6Ver = data.loc[:,'q6f']

    q1Ver = q1Ver[1:]
    q2Ver = q2Ver[1:]
    Sz = 0
    Sy = 1
    n1 =999
    n2= 1998
    n3=999
    n4 =1998
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4,x = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,Sz,Sy, n1,n2,n3,n4)

    q4VerM = np.concatenate((-q1Ver[::-1],-q6Ver[::-1]))
    q1VerM = (-q3Ver[::-1])
    q3VerM = (-q4Ver[::-1])
    q2VerM = np.concatenate((-q2Ver[::-1],-q5Ver))

    relError1 = ((q1[2,:]-q1VerM)/q1VerM) * 100
    relError2 = ((q2[2,:]-q2VerM)/q2VerM) * 100
    relError3 = ((q3[2,:]-q3VerM)/q3VerM) * 100
    relError4 = ((q4[2,:]-q4VerM)/q4VerM) * 100

    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Relative error VS distance in s-direction (Sy=1,Sz=0)", fontsize=16)

    axs[0, 0].plot(sVec1, relError1)
    axs[0, 0].set_title('Region 1')
    axs[0,0].grid()

    axs[0, 1].plot(sVec2, relError2)
    axs[0, 1].set_title('Region 2')
    axs[0,1].grid()

    axs[1, 0].plot(sVec3, relError3)
    axs[1, 0].set_title('Region 3')
    axs[1,0].grid()

    axs[1, 1].plot(sVec4, relError4)
    axs[1, 1].set_title('Region 4')
    axs[1,1].grid()
    axs[0,0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,0].set_xlabel('s [m]',fontsize = 12.0)
    axs[1, 0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,1].set_xlabel('s [m]',fontsize = 12.0)

    axs[0,0].tick_params(axis='both', which='major', labelsize=12)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)


    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Shear flow VS distance in s-direction (Sy=1,Sz=0)", fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    axs[0,0].plot(sVec1,q1[2,:],label="Numerical")
    axs[0,0].plot(sVec1,q1VerM,label="Verification")
    axs[0,0].set_title('Region 1')
    axs[0, 0].grid()

    axs[0, 1].plot(sVec2, q2[2,:], label="Numerical")
    axs[0, 1].plot(sVec2, q2VerM, label="Verification")
    axs[0, 1].set_title('Region 2')
    axs[0, 1].grid()

    axs[1, 0].plot(sVec3, q3[2,:], label="Numerical")
    axs[1, 0].plot(sVec3, q3VerM, label="Verification")
    axs[1, 0].set_title('Region 3')
    axs[1, 0].grid()

    axs[1, 1].plot(sVec4, q4[2,:], label="Numerical")
    axs[1, 1].plot(sVec4, q4VerM, label="Verification")
    axs[1, 1].set_title('Region 4')
    axs[1, 1].grid()

    axs[0, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    axs[1, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 1].set_xlabel('s [m]', fontsize=12.0)

    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.legend()
    plt.show()
    return

def shearFlowGraph2(craft):
    # draw graph for Sy = 0, Sz = 1 to compare with verification model
    #verification model shear flow, Sy = 0, Sz = 1
    data = pd.read_csv("verification_data/SzOne.csv", index_col=0)
    q1Ver = data.loc[:, 'q1f']
    q2Ver = data.loc[:, 'q2f']
    q3Ver = data.loc[:, 'q3f']
    q4Ver = data.loc[:, 'q4f']
    q5Ver = data.loc[:, 'q5f']
    q6Ver = data.loc[:, 'q6f']

    q1Ver = q1Ver[1:]
    q2Ver = q2Ver[1:]

    Sz = 1
    Sy = 0
    n1 =999
    n2=1998
    n3=999
    n4 =1998
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4,x = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,Sz,Sy, n1,n2,n3,n4)

    q4VerM = np.concatenate((-q1Ver[::-1],-q6Ver[::-1]))
    q1VerM = (-q3Ver[::-1])
    q3VerM = (-q4Ver[::-1])
    q2VerM = np.concatenate((-q2Ver[::-1],q5Ver))

    relError1 = np.divide((q1[2,:]-q1VerM), q1VerM, out=np.zeros_like((q1[2,:]-q1VerM)), where=q1VerM != 0) * 100
    relError2 = np.divide((q2[2,:] - q2VerM), q2VerM, out=np.zeros_like((q2[2,:] - q2VerM)), where=q2VerM != 0) * 100
    relError3 = np.divide((q3[2,:] - q3VerM), q3VerM, out=np.zeros_like((q3[2,:] - q3VerM)), where=q3VerM != 0) * 100
    relError4 = np.divide((q4[2,:] - q4VerM), q4VerM, out=np.zeros_like((q4[2,:] - q4VerM)), where=q4VerM != 0) * 100
    #relError1 = ((q1-q1VerM)/q1VerM) * 100
    #relError2 = ((q2-q2VerM)/q2VerM) * 100
    #relError3 = ((q3-q3VerM)/q3VerM) * 100
    #relError4 = ((q4-q4VerM)/q4VerM) * 100

    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Relative error VS distance in s-direction (Sy=0,Sz=1)", fontsize=16)

    axs[0, 0].plot(sVec1, relError1)
    axs[0, 0].set_title('Region 1')
    axs[0,0].grid()
    #axs[0,0].set_ylim(-5,5)

    axs[0, 1].plot(sVec2, relError2)
    axs[0, 1].set_title('Region 2')
    axs[0,1].grid()
    #axs[0, 1].set_ylim(-5, 5)

    axs[1, 0].plot(sVec3, relError3)
    axs[1, 0].set_title('Region 3')
    axs[1,0].grid()
    #axs[1, 0].set_ylim(-5, 5)

    axs[1, 1].plot(sVec4, relError4)
    axs[1, 1].set_title('Region 4')
    axs[1,1].grid()
    axs[0,0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,0].set_xlabel('s [m]',fontsize = 12.0)
    axs[1, 0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,1].set_xlabel('s [m]',fontsize = 12.0)

    axs[0,0].tick_params(axis='both', which='major', labelsize=12)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.show()


    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Shear flow VS distance in s-direction (Sy=0,Sz=1)", fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    axs[0,0].plot(sVec1,q1[2,:],label="Numerical")
    axs[0,0].plot(sVec1,q1VerM,label="Verification")
    axs[0,0].set_title('Region 1')
    axs[0, 0].grid()

    axs[0, 1].plot(sVec2, q2[2,:], label="Numerical")
    axs[0, 1].plot(sVec2, q2VerM, label="Verification")
    axs[0, 1].set_title('Region 2')
    axs[0, 1].grid()

    axs[1, 0].plot(sVec3, q3[2,:], label="Numerical")
    axs[1, 0].plot(sVec3, q3VerM, label="Verification")
    axs[1, 0].set_title('Region 3')
    axs[1, 0].grid()

    axs[1, 1].plot(sVec4, q4[2,:], label="Numerical")
    axs[1, 1].plot(sVec4, q4VerM, label="Verification")
    axs[1, 1].set_title('Region 4')
    axs[1, 1].grid()

    axs[0, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    axs[1, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 1].set_xlabel('s [m]', fontsize=12.0)

    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.legend()
    plt.show()
    return

def SySzHigh(craft):
    #draw graph for Sy = 1, Sz = 1 to compare with verification model
    #verification model shear flow, Sy = 1, Sz = 1
    data = pd.read_csv("verification_data/SzSyOne.csv", index_col=0)
    dataSz = pd.read_csv("verification_data/SzOne.csv", index_col=0)
    dataSy = pd.read_csv("verification_data/SyOne.csv", index_col=0)
    q1Ver = data.loc[:, 'q1f']
    q2Ver = data.loc[:, 'q2f']
    q3Ver = data.loc[:, 'q3f']
    q4Ver = data.loc[:, 'q4f']
    q6Ver = data.loc[:, 'q6f']

    q5Sz = dataSz.loc[:,'q5f']
    q5Sy = dataSy.loc[:, 'q5f']

    q1Ver = q1Ver[1:]
    q2Ver = q2Ver[1:]
    q5Ver = q5Sy - q5Sz

    Sz = 1
    Sy = 1
    n1 =999
    n2=1998
    n3=999
    n4 =1998
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4,x = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,Sz,Sy, n1,n2,n3,n4)



    q4VerM = np.concatenate((-q1Ver[::-1],-q6Ver[::-1]))
    q1VerM = (-q3Ver[::-1])
    q3VerM = (-q4Ver[::-1])
    q2VerM = np.concatenate((-q2Ver[::-1],-q5Ver))

    relError1 = np.divide((q1[2,:]-q1VerM), q1VerM, out=np.zeros_like((q1[2,:]-q1VerM)), where=q1VerM != 0) * 100
    relError2 = np.divide((q2[2,:] - q2VerM), q2VerM, out=np.zeros_like((q2[2,:] - q2VerM)), where=q2VerM != 0) * 100
    relError3 = np.divide((q3[2,:] - q3VerM), q3VerM, out=np.zeros_like((q3[2,:] - q3VerM)), where=q3VerM != 0) * 100
    relError4 = np.divide((q4[2,:] - q4VerM), q4VerM, out=np.zeros_like((q4[2,:] - q4VerM)), where=q4VerM != 0) * 100
    #relError1 = ((q1-q1VerM)/q1VerM) * 100
    #relError2 = ((q2-q2VerM)/q2VerM) * 100
    #relError3 = ((q3-q3VerM)/q3VerM) * 100
    #relError4 = ((q4-q4VerM)/q4VerM) * 100
    absError1 = np.average(np.abs(q1[2, :] - q1VerM))
    absError2 = np.average(np.abs(q2[2, :] - q2VerM))
    absError3 = np.average(np.abs(q3[2, :] - q3VerM))
    absError4 = np.average(np.abs(q4[2, :] - q4VerM))
    print("absError1 high", absError1)
    print("absError2 high", absError2)
    print("absError3 high", absError3)
    print("absError4 high", absError4)

    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Relative error VS distance in s-direction (Sy=1,Sz=1) (n1=%d,n2=%d,n3=%d,n4=%d)"%(n1+1,n2+1,n3+1,n4+1), fontsize=16)

    axs[0, 0].plot(sVec1, relError1)
    axs[0, 0].set_title('Region 1')
    axs[0,0].grid()
    axs[0,0].set_ylim(-1,1)

    axs[0, 1].plot(sVec2, relError2)
    axs[0, 1].set_title('Region 2')
    axs[0,1].grid()
    axs[0, 1].set_ylim(0, 0.5)

    axs[1, 0].plot(sVec3, relError3)
    axs[1, 0].set_title('Region 3')
    axs[1,0].grid()
    axs[1, 0].set_ylim(-1, 1)

    axs[1, 1].plot(sVec4, relError4)
    axs[1, 1].set_title('Region 4')
    axs[1, 1].set_ylim(-1, 1)
    axs[1,1].grid()
    axs[0,0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,0].set_xlabel('s [m]',fontsize = 12.0)
    axs[1, 0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,1].set_xlabel('s [m]',fontsize = 12.0)

    axs[0,0].tick_params(axis='both', which='major', labelsize=12)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)


    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Shear flow VS distance in s-direction (Sy=1,Sz=1) (n1=%d,n2=%d,n3=%d,n4=%d)"%(n1+1,n2+1,n3+1,n4+1), fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    axs[0,0].plot(sVec1,q1[2,:],label="Numerical")
    axs[0,0].plot(sVec1,q1VerM,label="Verification")
    axs[0,0].set_title('Region 1')
    axs[0, 0].grid()

    axs[0, 1].plot(sVec2, q2[2,:], label="Numerical")
    axs[0, 1].plot(sVec2, q2VerM, label="Verification")
    axs[0, 1].set_title('Region 2')
    axs[0, 1].grid()

    axs[1, 0].plot(sVec3, q3[2,:], label="Numerical")
    axs[1, 0].plot(sVec3, q3VerM, label="Verification")
    axs[1, 0].set_title('Region 3')
    axs[1, 0].grid()

    axs[1, 1].plot(sVec4, q4[2,:], label="Numerical")
    axs[1, 1].plot(sVec4, q4VerM, label="Verification")
    axs[1, 1].set_title('Region 4')
    axs[1, 1].grid()

    axs[0, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    axs[1, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 1].set_xlabel('s [m]', fontsize=12.0)

    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)

    #fig, axs = plt.subplots(2, 2)
    #plt.suptitle("Absolute error VS distance in s-direction (Sy=1,Sz=1) (n1=%d,n2=%d,n3=%d,n4=%d)" % (
    #n1 + 1, n2 + 1, n3 + 1, n4 + 1), fontsize=16)
    #plt.rcParams["font.size"] = "12"
    #plt.rcParams["axes.labelsize"] = "12"
    #axs[0, 0].plot(sVec1, absError1)
    #axs[0, 0].set_title('Region 1')
    #axs[0, 0].grid()
#
    #axs[0, 1].plot(sVec2, absError2)
    #axs[0, 1].set_title('Region 2')
    #axs[0, 1].grid()
#
    #axs[1, 0].plot(sVec3, absError3)
    #axs[1, 0].set_title('Region 3')
    #axs[1, 0].grid()
#
    #axs[1, 1].plot(sVec4, absError4)
    #axs[1, 1].set_title('Region 4')
    #axs[1, 1].grid()
#
    #axs[0, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    #axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    #axs[1, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    #axs[1, 1].set_xlabel('s [m]', fontsize=12.0)
#
    #axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    #axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    #axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    #axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    #axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    #axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    #axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    #axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)


    plt.legend()
    #plt.show()
    return

def SySzLow(craft):
    #draw graph for Sy = 1, Sz = 1 to compare with verification model
    #verification model shear flow, Sy = 1, Sz = 1
    data = pd.read_csv("verification_data/SzSyOneLOW.csv", index_col=0)
    dataSz = pd.read_csv("verification_data/SzLOW.csv", index_col=0)
    dataSy = pd.read_csv("verification_data/SyLOW.csv", index_col=0)
    q1Ver = data.loc[:, 'q1f']
    q2Ver = data.loc[:, 'q2f']
    q3Ver = data.loc[:, 'q3f']
    q4Ver = data.loc[:, 'q4f']
    q6Ver = data.loc[:, 'q6f']

    q5Sz = dataSz.loc[:,'q5f']
    q5Sy = dataSy.loc[:, 'q5f']

    q1Ver = q1Ver[1:]
    q2Ver = q2Ver[1:]
    q5Ver = q5Sy - q5Sz

    Sz = 1
    Sy = 1
    n1 =99
    n2=198
    n3=99
    n4 =198
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4,x = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,Sz,Sy, n1,n2,n3,n4)



    q4VerM = np.concatenate((-q1Ver[::-1],-q6Ver[::-1]))
    q1VerM = (-q3Ver[::-1])
    q3VerM = (-q4Ver[::-1])
    q2VerM = np.concatenate((-q2Ver[::-1],-q5Ver))

    relError1 = np.divide((q1[2,:]-q1VerM), q1VerM, out=np.zeros_like((q1[2,:]-q1VerM)), where=q1VerM != 0) * 100
    relError2 = np.divide((q2[2,:] - q2VerM), q2VerM, out=np.zeros_like((q2[2,:] - q2VerM)), where=q2VerM != 0) * 100
    relError3 = np.divide((q3[2,:] - q3VerM), q3VerM, out=np.zeros_like((q3[2,:] - q3VerM)), where=q3VerM != 0) * 100
    relError4 = np.divide((q4[2,:] - q4VerM), q4VerM, out=np.zeros_like((q4[2,:] - q4VerM)), where=q4VerM != 0) * 100
    #relError1 = ((q1-q1VerM)/q1VerM) * 100
    #relError2 = ((q2-q2VerM)/q2VerM) * 100
    #relError3 = ((q3-q3VerM)/q3VerM) * 100
    #relError4 = ((q4-q4VerM)/q4VerM) * 100
    absError1 = np.average(np.abs(q1[2, :] - q1VerM))
    absError2 = np.average(np.abs(q2[2, :] - q2VerM))
    absError3 = np.average(np.abs(q3[2, :] - q3VerM))
    absError4 = np.average(np.abs(q4[2, :] - q4VerM))
    print("absError1 low",absError1)
    print("absError2 low", absError2)
    print("absError3 low", absError3)
    print("absError4 low", absError4)

    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Relative error VS distance in s-direction (Sy=1,Sz=1) (n1=%d,n2=%d,n3=%d,n4=%d)"%(n1+1,n2+1,n3+1,n4+1), fontsize=16)

    axs[0, 0].plot(sVec1, relError1)
    axs[0, 0].set_title('Region 1')
    axs[0,0].grid()
    axs[0,0].set_ylim(-1,1)

    axs[0, 1].plot(sVec2, relError2)
    axs[0, 1].set_title('Region 2')
    axs[0,1].grid()
    axs[0, 1].set_ylim(0, 0.5)

    axs[1, 0].plot(sVec3, relError3)
    axs[1, 0].set_title('Region 3')
    axs[1,0].grid()
    axs[1, 0].set_ylim(-1, 1)

    axs[1, 1].plot(sVec4, relError4)
    axs[1, 1].set_title('Region 4')
    axs[1, 1].set_ylim(-1, 1)
    axs[1,1].grid()
    axs[0,0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,0].set_xlabel('s [m]',fontsize = 12.0)
    axs[1, 0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,1].set_xlabel('s [m]',fontsize = 12.0)

    axs[0,0].tick_params(axis='both', which='major', labelsize=12)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)


    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Shear flow VS distance in s-direction (Sy=1,Sz=1) (n1=%d,n2=%d,n3=%d,n4=%d)"%(n1+1,n2+1,n3+1,n4+1), fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    axs[0,0].plot(sVec1,q1[2,:],label="Numerical")
    axs[0,0].plot(sVec1,q1VerM,label="Verification")
    axs[0,0].set_title('Region 1')
    axs[0, 0].grid()

    axs[0, 1].plot(sVec2, q2[2,:], label="Numerical")
    axs[0, 1].plot(sVec2, q2VerM, label="Verification")
    axs[0, 1].set_title('Region 2')
    axs[0, 1].grid()

    axs[1, 0].plot(sVec3, q3[2,:], label="Numerical")
    axs[1, 0].plot(sVec3, q3VerM, label="Verification")
    axs[1, 0].set_title('Region 3')
    axs[1, 0].grid()

    axs[1, 1].plot(sVec4, q4[2,:], label="Numerical")
    axs[1, 1].plot(sVec4, q4VerM, label="Verification")
    axs[1, 1].set_title('Region 4')
    axs[1, 1].grid()

    axs[0, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    axs[1, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    axs[1, 1].set_xlabel('s [m]', fontsize=12.0)

    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)

    #fig, axs = plt.subplots(2, 2)
    #plt.suptitle("Absolute error VS distance in s-direction (Sy=1,Sz=1) (n1=%d,n2=%d,n3=%d,n4=%d)" % (
    #    n1 + 1, n2 + 1, n3 + 1, n4 + 1), fontsize=16)
    #plt.rcParams["font.size"] = "12"
    #plt.rcParams["axes.labelsize"] = "12"
    #axs[0, 0].plot(sVec1, absError1)
    #axs[0, 0].set_title('Region 1')
    #axs[0, 0].grid()
#
    #axs[0, 1].plot(sVec2, absError2)
    #axs[0, 1].set_title('Region 2')
    #axs[0, 1].grid()
#
    #axs[1, 0].plot(sVec3, absError3)
    #axs[1, 0].set_title('Region 3')
    #axs[1, 0].grid()
#
    #axs[1, 1].plot(sVec4, absError4)
    #axs[1, 1].set_title('Region 4')
    #axs[1, 1].grid()
#
    #axs[0, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    #axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    #axs[1, 0].set_ylabel('Shear flow [N\m]', fontsize=12.0)
    #axs[1, 1].set_xlabel('s [m]', fontsize=12.0)
#
    #axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    #axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    #axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    #axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    #axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    #axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    #axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    #axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.legend()
    #plt.show()
    return