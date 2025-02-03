# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 16:02:04 2024
@author: mendez
"""

''' 
This file generates the data for the Pendulum assimilation problem
We simulate some noisy measurements on theta and theta dot in the pendulum
exercise.

'''

import numpy as np
import matplotlib.pyplot as plt
# We use odeint for the time integration
from scipy.integrate import odeint
from scipy.stats import uniform # to generate uniform distr


# Customize plot settings for LaTeX and larger fonts
plt.rcParams.update({
    'text.usetex': True,
    'font.size': 16,
    'font.family': 'serif'
})


# We define the function of the system    
def f(s, t,mu=0.5,omega_n=5):
    # this is the function f in the slides
    s_1_dot = s[1]
    s_2_dot = -mu*s[1]-omega_n**2*np.sin(s[0])
    return [s_1_dot, s_2_dot]

# To make odeint go faster, we provide the Jacobian    
def df_ds(s,t,mu,omega_n):
    s_1,s_2=s
    jac=[[0,1],[-omega_n**2*np.cos(s_1),-mu]]
    return jac




#%% Data about the pendulum (educated guess)

m=120/1000 # mass of the pendulum  [kg]
L=41/100 # full length of the pendulum  [m]
w=37/1000 # width of the pendulum [m]
l_cm=24/100 # distance from center of mass to pivot [m]
g=9.815 # gravitational acceleration [m/s^2]
I_x=1/12*m*(L**2+w**2)+m*l_cm**2 # Moment of Inertia [kg m^2]

omega_n=np.sqrt(m*g*l_cm/I_x) # natural omega
T_n=2*np.pi/omega_n # natural period [s]


degs=np.array([+35,-28,+15,-45,36,-22,18,42])
s0_s=degs*np.pi/180
mu_s=np.array([0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18])

import os
FOLDER_DATA='data_pendulum'
if not os.path.exists(FOLDER_DATA):
    os.mkdir(FOLDER_DATA)

for k in range(8):
    
    # Define the file name
    file_name = FOLDER_DATA+os.sep+'measurement_'+str(k)+'.dat'
    s0=[s0_s[k],0]
    
    mu=mu_s[k] # very much unknown!     


    n_t=501; 
    T_0,T_f=0,10
    t=np.linspace(T_0,T_f,n_t); 
    
    Y = odeint(f, s0, t,args=(mu,omega_n))
    s1_f=Y[:,0]; s2_f=Y[:,1] # Get the solution
    
    
    # Add noise and spikes
    noise=s1_f.max()/3*np.random.rand(n_t)
    spikes=np.zeros(n_t);
    
    for j in range(n_t//20):
        # pick a random index in 0,n_t-1
        INDEX=np.random.randint(0, n_t-1)
        # pick a random continuos value between -0.4 and 0.4
        SPIKE=uniform.rvs(size=1, loc = -s1_f.max()/2, scale=s1_f.max())
        spikes[INDEX]=SPIKE
    
    
    s1_m=s1_f+noise+spikes
    
    # Here's the plot of the solution    
    fig, ax = plt.subplots(figsize=(6,4)) # Create Signal Noisy 
    plt.plot(t,s1_m)
    plt.xlabel('$t$',fontsize=18)
    plt.ylabel('$\\theta(t)$',fontsize=18)
    #plt.title('Regression')
    plt.tight_layout()
    plt.show()
    Name=FOLDER_DATA+os.sep+'Theta_Evol_'+str(k)+'.png'
    plt.savefig(Name, dpi=300)      
    # plt.close(fig)
    
    
    # Add noise and spikes
    noise=s2_f.max()/3*np.random.rand(n_t)
    spikes=np.zeros(n_t);
    
    for j in range(n_t//20):
        # pick a random index in 0,n_t-1
        INDEX=np.random.randint(0, n_t-1)
        # pick a random continuos value between -0.4 and 0.4
        SPIKE=uniform.rvs(size=1, loc = -s2_f.max(), scale=s2_f.max())
        spikes[INDEX]=SPIKE
    
    
    s2_m=s2_f+noise+spikes
    
    # Here's the plot of the solution    
    fig, ax = plt.subplots(figsize=(6,4)) # Create Signal Noisy 
    plt.plot(t,s2_m)
    plt.xlabel('$t$',fontsize=18)
    plt.ylabel('$\dot{\\theta}(t)$',fontsize=18)
    # plt.legend(shadow=True,fontsize=17)
    #plt.title('Regression')
    plt.tight_layout()
    plt.show()
    Name=FOLDER_DATA+os.sep+'Theta_dot_Evol_'+str(k)+'.png'
    plt.savefig(Name, dpi=300)      
    # plt.close(fig)
    
    
    # Export the data in a dat file 
    
    
    # Combine the arrays into a single 2D array
    data = np.column_stack((t.astype(np.float16), s1_m.astype(np.float16), s2_m.astype(np.float16)))
    
    
    # Save the data to a .dat file
    np.savetxt(file_name, data, header='time theta theta_dot', comments='', delimiter='\t',fmt='%1.7f')
    
    print(f'Data saved to {file_name}')

























