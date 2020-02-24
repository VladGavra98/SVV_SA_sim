# -*- coding: utf-8 -*-
"""
Calculates the normal stresses and outputs an array of maximum stresses and their location for each section

@author: alberto
"""

def calcNormStress(Ha,Ca,Nst,MomentArray,StPos,Zcg,Iyy,Izz):
    # Only 2nd and 3rd row of MomentArray are used
    
    
    # Append spar cap locations to StPos
    CapsPos = np.array([[Ha/2,-Ha/2],
                        [-Ha/2,-Ha/2]])
    LocationArray = np.append(StPos,CapsPos,axis=1)
    
    Nsections = len(MomentArray[0,:]) #amount of cross sections
    Npoints = len(LocationArray[0,:]) #amount of points
    
    # ---------------- Sigma -------------------------
    # the goal is to create an array with only the maximum stress and its location for each section
    # This array is MaxStresses
    MaxStresses = np.zeros((3,Nsections))   #array of shape (y,z,sigma_max) per section
    LocalStresses = np.zeros((3,Npoints))   #array of shape (y,z,sigma_max) per point
    
    for i in range(Nsections):                                          # i covers each cross section
        for n in range(Npoints):                                        # n covers each point
            sigma_y = MomentArray[1,i]*(LocationArray[1,n]-Zcg)/Iyy     # compute stress
            sigma_z = -MomentArray[2,i]*LocationArray[0,n]/Izz
            LocalStresses[0,n] = LocationArray[0,n]                     # add to local list
            LocalStresses[1,n] = LocationArray[1,n]
            LocalStresses[2,n] = sigma_y + sigma_z
              
        #print(LocalStresses)
        LocalMax = np.amax(np.abs(LocalStresses), axis=1)[2]        # find maximum stress (= third row, index 2)
        #print("Local Max", LocalMax)
        IndexMax = np.argmax(np.abs(LocalStresses),axis=1)[2]     # find its location
        #print("Index Max:",IndexMax)
        MaxStresses[0,i] = LocationArray[0,IndexMax]                  # store its y coordinate
        MaxStresses[1,i] = LocationArray[1,IndexMax]                  # store its z coordinate
        MaxStresses[2,i] = LocalStresses[2,IndexMax]
    
    return MaxStresses