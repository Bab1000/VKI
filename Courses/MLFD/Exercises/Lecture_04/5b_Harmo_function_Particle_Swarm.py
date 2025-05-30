# -*- coding: utf-8 -*-
"""
Created on Jan 23 15:03:41 2022

@author: mendez
"""


import numpy as np
import matplotlib.pyplot as plt


#%% Preamble: customization of matplotlib
# Configuration for plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

#%% Harmonic function definition
def harmo(X):  # The Function to Minimize  
    """Function to be minimized"""  # new        
    C=X[0]*np.sin(4*X[0])+1.1*X[1]*np.sin(2*X[1])
    return C 


from PSO_Functions import Initialize_POP
from PSO_Functions import Evaluate_POP
from PSO_Functions import Update_POP
from PSO_Functions import PSO
from PSO_Functions import Anim_COMP

#---------------- Parameters
N_POP=400
X_Bounds=[(-8,8),(-8,8)]
Func=harmo
N_ITER=50


X_V, V_P=Initialize_POP(N_POP,X_Bounds,n_G=0.5,sigma_I_r=6,I_V=0.1)
Err_1=Evaluate_POP(X_V,Func)

X_V_n, V_P_n, X_B_V, Err_B_V=Update_POP(X_V,V_P,X_V,Err_1,Err_1,
               X_Bounds,n_I=0,N_ITER=100,w_I=0.3,w_F=0.05,c_c=2,c_s=2)

# Example of non collaborative particles
# X_S, X_U, X_V=PSO(Func,X_Bounds,n_p=400,N_ITER=400,
#        n_G=0.5,sigma_I_r=6,w_I=0.0001,w_F=0.0001,c_c=0.1,c_s=0.3)

# Example with collaborative particles
X_S, X_U, X_V=Anim_COMP(Func,X_Bounds,n_p=100,N_ITER=100,\
                 n_G=0.5,sigma_I_r=6,w_I=0.3,w_F=0.01,c_c=0.02,c_s=0.5,\
                 x_1m=-8,x_1M=8,x_2m=-8,x_2M=8,\
                     npoints=200,Name_Video='Harmo_PSO_TW.gif')

    
# Example with non collaborative particles
X_S, X_U, X_V=Anim_COMP(Func,X_Bounds,n_p=100,N_ITER=100,\
                 n_G=0.5,sigma_I_r=6,w_I=0.3,w_F=0.01,c_c=0.8,c_s=0.02,\
                 x_1m=-8,x_1M=8,x_2m=-8,x_2M=8,\
                     npoints=200,Name_Video='Harmo_PSO_EG.gif')
    