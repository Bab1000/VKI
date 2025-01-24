

import numpy as np
import csv
import math
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def fit_parabola(x, y):
    a, b, c = np.polyfit(x, y, 2)
    max_x = -b/(2*a)
    return max_x

x_probe=0.871                     # Probe location [m]
x_Torchexit=0.486
r_probe=0.025                     # Probe radius [m]
mdot=16                           # Mass flow rate [g/s]  
ps= 50 * 100                     # Static pressure [mbar --> Pa]

# Powers kW
P=[100]
#P=[125, 130, 140, 145, 170, 180, 190]
#P=[50, 55, 60, 75, 80, 90, 100, 110, 120, 125, 130, 140, 145, 150, 160, 170, 175, 180, 190, 200]     # 10 mbar
#P=[50, 55, 60, 70, 75, 80, 90, 100, 110, 120, 130, 140, 145, 150, 157, 165, 176, 187, 198, 200]      # 30 mbar
#P=[50, 55, 60, 65, 70, 75, 78, 80, 85, 88, 90, 92, 95, 98, 100, 102, 105, 108, 110, 113, 115, 117, 120, 125, 127, 130, 133, 135, 140, 145, 148, 150, 155, 160, 165, 170, 173, 175, 180, 185, 190, 195, 200]   # 50 mbar
#P=[50, 55, 60, 70, 75, 80, 90, 100, 110, 120, 130, 141, 145, 150, 160, 170, 180, 190, 200]           #75 mbar
#P=[50, 55, 60, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 132, 135, 140, 143, 145, 150, 155, 160, 170, 180, 190, 200]        # 100 mbar


# Path to files
#os.chdir('/mnt/c/Users/xavie/Documents/Unif Justin/VKI/Project/ICP_Roemer/results/full_sims')
path_to_data = '/mnt/c/Users/xavie/Documents/Unif Justin/VKI/Project/ICP_Roemer/results/full_sims'

# Create csv file

with open('Results/FS_15mbar_80kW_16gs_air.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=';')
    writer.writerow(['mdot [g/s]', 'ps [mbar]', 'Pel [kW]', 'u_Torchexit [m/s]', 'beta_Torchexit [1/s]', 
                     'T_Torchexit [K]', 'p_Torchexit [Pa]', 'rho_Torchexit [kg/m^3]', 'h_Torchexit [J/kg]', 
                     'u_freestream [m/s]', 'beta_freestream [1/s]', 'T_freestream [K]', 'p_freestream [Pa]', 
                     'rho_freestream [kg/m^3]', 'h_freestream [J/kg]', 'delta [m]', 'u_delta [m/s]', 
                     'beta_delta [1/s]', 'dbetadx_delta [1/ms]', 'T_delta [K]', 'p_delta [Pa]', 
                     'rho_delta [kg/m^3]', 'h_delta [J/kg]', 'u_wall [m/s]', 'beta_wall [1/s]', 
                     'T_wall [K]', 'p_wall [Pa]', 'rho_wall [kg/m^3]', 'h_wall [J/kg]', 
                     'NDP1', 'NDP2', 'NDP3', 'NDP4', 'NDP5'])


for i in P:
    Pel=i                            # Power [kW]

    name_folder= path_to_data + "/FS_15mbar_80kW_16gs_air.plt"     # Name of the folder of simulation

    WithProbe = pd.read_table(name_folder, sep='\s+', skiprows=np.arange(3),  names=['x', 'y', 'dp', 'u', 'v', 'Te', 'EpR', 'EpI', 'rho', 'H', 'M', 'xc0', 'xc1', 'xc2', 'xc3', 'xc4', 'xc5', 'xc6', 'xc7', 'xc8', 'xc9', 'xc10'])
    x=np.array(WithProbe.x)
    y=np.array(WithProbe.y)  
    u=np.array(WithProbe.u)
    v=np.array(WithProbe.v)
    dp=np.array(WithProbe.dp)
    T=np.array(WithProbe.Te)
    H=np.array(WithProbe.H)
    rho=np.array(WithProbe.rho)
    M=np.array(WithProbe.M)
    N=len(x)
    x_axis=[]
    u_axis=[]
    v_axis=[]
    dp_axis=[]
    T_axis=[]
    H_axis=[]
    rho_axis=[]
    M_axis=[]
    u_axis_sorted=[]
    v_axis_sorted=[]
    dp_axis_sorted=[]
    T_axis_sorted=[]
    H_axis_sorted=[]
    rho_axis_sorted=[]
    M_axis_sorted=[]
    for i in range(N):
        if y[i]==0:
            x_axis.append(x[i])
            u_axis.append(u[i])
            v_axis.append(v[i])
            dp_axis.append(dp[i])
            T_axis.append(T[i])
            H_axis.append(H[i])
            rho_axis.append(rho[i])
            M_axis.append(M[i])
        
    x_axis_sorted=sorted(x_axis)
    NN=len(x_axis_sorted)
    for i in range(NN):
        for j in range(NN):
            if x_axis[j]==x_axis_sorted[i]:
                u_axis_sorted.append(u_axis[j])
                v_axis_sorted.append(v_axis[j])
                dp_axis_sorted.append(dp_axis[j])
                T_axis_sorted.append(T_axis[j])
                H_axis_sorted.append(H_axis[j])
                rho_axis_sorted.append(rho_axis[j])
                M_axis_sorted.append(M_axis[j])

################################### Assessement of the velocity gradient 'Beta' ###############################################################

# First find the first points with y-coordinate different from zero
    Torchpoints = np.count_nonzero(np.array(x_axis_sorted) <= 0.486)                                                                                # Number of points inside the torch
    bool_arr = np.logical_and(np.array(x_axis_sorted)>=0.851, np.array(x_axis_sorted) <= 0.871)
    BLpoints = np.count_nonzero(bool_arr)                                                                                                           # Number of points inside the mesh boundary layer 
    y1Torch=0.002                                                                                                                                   # Distance in meter of the first node from the axis in the torch
    y1Probe=0.0001                                                                                                                                  # Distance in meter of the first node from the axis in the mesh boundary layer
    y1TorchProbe=np.interp(x_axis_sorted[Torchpoints:len(x_axis_sorted)-BLpoints], [0.486, 0.851], [y1Torch, y1Probe], left=None, right=None)       # Estimated distance in meter of the first node from the axis in the zone between the Torch and mesh boundary layer
    y1prior=np.concatenate((np.ones(Torchpoints)*y1Torch, y1TorchProbe))
    y1prior=np.concatenate((y1prior,np.ones(BLpoints)*y1Probe))                                                                                     # Estimated distance in meter of the first node from the axis along the stagnation line
    closest_point_index=[]
    for j in range(NN):
        target_coords = [x_axis_sorted[j], y1prior[j]]
        distance=[]
        for i in range(N):
            distance.append(np.sqrt((x[i]-target_coords[0])**2+(y[i]-target_coords[1])**2))
        closest_point_index.append(np.argmin(distance))

# Coordinate of the first point with y-coordinate different from zero 
    x1=[]
    y1=[]
    v1=[]
    for i in range(NN):
        x1.append(x[closest_point_index[i]])
        y1.append(y[closest_point_index[i]])
        v1.append(v[closest_point_index[i]])
    v1_interp=np.interp(x_axis_sorted, x1, v1, left=None, right=None)
    y1_interp=np.interp(x_axis_sorted, x1, y1, left=None, right=None)

# Assessmente of Beta and its derivative with respect to x

    beta=v1_interp/(np.array(y1_interp)+1e-6)
    dbetadx=np.zeros(NN)
    for i in np.arange(1, NN-1):
        dbetadx[i] = (beta[i+1]-beta[i-1])/(x_axis_sorted[i+1]-x_axis_sorted[i-1])                           # Second order first central derivative
    dbetadx[0] = (-3*beta[0]+4*beta[1]-beta[2])/(x_axis_sorted[1]-x_axis_sorted[0])                         # Second order first forward derivative
    dbetadx[NN-1] = (3*beta[NN-1]-4*beta[NN-2]+beta[NN-3])/(x_axis_sorted[NN-1]-x_axis_sorted[NN-2])                  # Second order first backword derivative
#    max_indices = argrelextrema(dbetadx, np.greater)
#    BL_index=max_indices[-1][-1]
    max_indices = argrelextrema(dbetadx, np.greater)
    n_maxindices = len(max_indices[0])
    BL_index=max_indices[0][n_maxindices-2]
# Assessment of the exact x-value where there is the inflection point of beta
# I fit a parabola on three indices (that where there is the relative maximum plus the previous and the succesive one)

    x_fit = x_axis_sorted[BL_index-1:BL_index+2]                # x-values on which I will fit the parabola
    y_fit = dbetadx[BL_index-1:BL_index+2]                      # Parabola around the index of relative maximum of beta
    x_BLfitted = fit_parabola(x_fit, y_fit)                     # Exact x-value of the boundary layer edge

################################ ICP input ################################################

#print("\n ICP input: \n mdot = {} g/s \n ps = {} mbar \n Pel = {} kW \n".format(mdot, math.trunc(ps/100), Pel))

################################ Torch exit ###############################################

    u_Torchexit=np.interp(x_Torchexit, x_axis_sorted, u_axis_sorted, left=None, right=None)
    beta_Torchexit=np.interp(x_Torchexit, x_axis_sorted, beta, left=None, right=None)
    T_Torchexit=np.interp(x_Torchexit, x_axis_sorted, T_axis_sorted, left=None, right=None)
    p_Torchexit=ps+np.interp(x_Torchexit, x_axis_sorted, dp_axis_sorted, left=None, right=None)
    H_Torchexit=np.interp(x_Torchexit, x_axis_sorted, H_axis_sorted, left=None, right=None)
    rho_Torchexit=np.interp(x_Torchexit, x_axis_sorted, rho_axis_sorted, left=None, right=None)
    M_Torchexit=np.interp(x_Torchexit, x_axis_sorted, M_axis_sorted, left=None, right=None)

#print("\n Torch exit: \n x = {} m \n u = {} m/s \n beta = {} 1/s \n T = {} K \n p = {} Pa \n rho = {} kg/m^3 \n h = {} J/kg \n".format(x_Torchexit, u_Torchexit, beta_Torchexit, T_Torchexit, p_Torchexit, rho_Torchexit, H_Torchexit))

############################### Freestream point ##########################################

    x_freestream=x_Torchexit+(x_probe-x_Torchexit)/2
    u_freestream=np.interp(x_freestream, x_axis_sorted, u_axis_sorted, left=None, right=None)
    beta_freestream=np.interp(x_freestream, x_axis_sorted, beta, left=None, right=None)
    T_freestream=np.interp(x_freestream, x_axis_sorted, T_axis_sorted, left=None, right=None)
    p_freestream=ps+np.interp(x_freestream, x_axis_sorted, dp_axis_sorted, left=None, right=None)
    H_freestream=np.interp(x_freestream, x_axis_sorted, H_axis_sorted, left=None, right=None)
    rho_freestream=np.interp(x_freestream, x_axis_sorted, rho_axis_sorted, left=None, right=None)
    M_freestream=np.interp(x_freestream, x_axis_sorted, M_axis_sorted, left=None, right=None)

#print("\n Freestream point: \n x = {} m \n u = {} m/s \n beta = {} 1/s \n T = {} K \n p = {} Pa \n rho = {} kg/m^3 \n h = {} J/kg \n".format(x_freestream, u_freestream, beta_freestream, T_freestream, p_freestream, rho_freestream, H_freestream))

################################ Boundary Layer Quantities ################################

    delta_BL=x_probe-x_BLfitted
    u_BL=np.interp(x_BLfitted, x_axis_sorted, u_axis_sorted, left=None, right=None)
    beta_BL=np.interp(x_BLfitted, x_axis_sorted, beta, left=None, right=None)
    dbetadx_BL=np.interp(x_BLfitted, x_axis_sorted, dbetadx, left=None, right=None)
    p_BL=ps+np.interp(x_BLfitted, x_axis_sorted, dp_axis_sorted, left=None, right=None)
    T_BL=np.interp(x_BLfitted, x_axis_sorted, T_axis_sorted, left=None, right=None)
    H_BL=np.interp(x_BLfitted, x_axis_sorted, H_axis_sorted, left=None, right=None)
    rho_BL=np.interp(x_BLfitted, x_axis_sorted, rho_axis_sorted, left=None, right=None)
    M_BL=np.interp(x_BLfitted, x_axis_sorted, M_axis_sorted, left=None, right=None)

#print("\n Boundary Layer edge: \n x = {} m \n delta = {} m \n u = {} m/s \n beta = {} 1/s \n dbeta_dx = {} 1/(m s) \n T = {} K \n p = {} Pa \n rho = {} kg/m^3 \n h = {} J/kg \n".format(x_BLfitted, delta_BL, u_BL, beta_BL, dbetadx_BL, T_BL, p_BL, rho_BL, H_BL))

############################## Wall Point ##################################################

    u_wall=u_axis_sorted[-1]
    beta_wall=beta[-1]
    T_wall=T_axis_sorted[-1]
    p_wall=ps+dp_axis_sorted[-1]
    H_wall=H_axis_sorted[-1]
    rho_wall=rho_axis_sorted[-1]
    M_wall=M_axis_sorted[-1]

#print("\n Wall point: \n x = {} m \n u = {} m/s \n beta = {} 1/s \n T = {} K \n p = {} Pa \n rho = {} kg/m^3 \n h = {} J/kg \n".format(x_probe, u_wall, beta_wall, T_wall, p_wall, rho_wall, H_wall))

############################## Non-Dimensional Parameters ##################################

    NDP1=delta_BL/r_probe
    NDP2=beta_BL*r_probe/u_Torchexit
    NDP3=r_probe**2/u_Torchexit*dbetadx_BL
    NDP4=u_BL/u_Torchexit
    NDP5=u_BL/u_freestream

#print("\n Non-Dimensional Parameters: \n NDP1 = {} \n NDP2 = {} \n NDP3 = {} \n NDP4 = {} \n NDP5 = {} \n".format(NDP1, NDP2, NDP3, NDP4, NDP5))

################################ Write in a file ###########################################   
    result=[mdot, math.trunc(ps/100), Pel, u_Torchexit, beta_Torchexit, T_Torchexit, p_Torchexit, rho_Torchexit, H_Torchexit, u_freestream, beta_freestream, T_freestream, p_freestream, rho_freestream, H_freestream, delta_BL, u_BL, beta_BL, dbetadx_BL, T_BL, p_BL, rho_BL, H_BL, u_wall, beta_wall, T_wall, p_wall, rho_wall, H_wall, NDP1, NDP2, NDP3, NDP4, NDP5]
    
    f = open('Results/FS_15mbar_80kW_16gs_air.csv', 'a')
    writer = csv.writer(f, delimiter=';')
    writer.writerow(result)
    f.close()
    
############################## Plot ########################################################

    
    plt.close()
    plt.rc('text', usetex=True)       # This is for plot customization
    plt.rc('font', family='serif')
    plt.rc('lines', linewidth = 3)
    plt.rc('xtick',labelsize=18)
    plt.rc('ytick',labelsize=18)
    plt.rc('axes', labelsize=20)      # fontsize of the x and y labels 


    plt.figure()
    plt.suptitle(r'$p_s$ = {} mbar $\|$ $P$ = {} kW '.format(math.trunc(ps/100), Pel), fontsize=20)
    plt.subplot(2,1,1)
    plt.plot(x_axis_sorted,beta/10**4)
    #plt.plot([x_BLfitted, x_BLfitted], [np.min(beta / 10**4), np.max(beta / 10**4)], '--r', linewidth=2)
    # plt.xlabel('x [m]')
    plt.ylabel(r'$\beta$ [$10^4$ 1/s]')
    plt.xlim(min(0.82, x_BLfitted), max(0.890, x_BLfitted))
    plt.ylim(0,np.max(beta/10**4)+0.5)
    # plt.grid()
    plt.xticks([])                                  # Hide y axis
    # plt.tick_params(axis='x', length=0)             #
    # plt.tight_layout()
    
    plt.subplot(2,1,2)
    plt.plot(x_axis_sorted, dbetadx/10**6)
    plt.plot([x_BLfitted, x_BLfitted], [np.min(dbetadx/10**6), np.max(dbetadx/10**6)], '--k')
    plt.xlabel('x [m]')
    plt.ylabel(r'd$\beta$/dx [$10^6$ 1/(m$\cdot$s)]')
    plt.xlim(0.82,0.890)
    plt.ylim(-1, dbetadx[BL_index]/10**6+      1)   # changed 0.5 to 1 
    # plt.grid()
    plt.tight_layout()
    #plt.show()
    plt.savefig('Results/Beta/Boundarylayeredge_' + str(math.trunc(ps/100)) + 'mbar_' + str(Pel) + 'kW.png')
    plt.close('all')