# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:02:51 2023

@author: mendez
"""

# We solve the same problem in files 1 and 2 using a backward integraion.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy import optimize
import time


# Configuration for plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

# Here is the definition of the system with the default parameters
a,b,c,d=1,1,3,1

# Define the 'Good parameters':
p_star=np.array([a,b,c,d])

# Define the 'Initial Guesses'
p0=np.array([1.2,1.4,3.4,1.2])

def f(t, u, a=1, b=1, c=3, d=1):
    u1, u2 = u
    return [a*u1 - b*u1*u2, -c*u2 + d*u1*u2]

# Load the data
data = np.load('Data_Ex_Assimilation.npz')
t_d=data['t']
u_d=data['u_d']

#%% Small upgrade :) ! we use an interpolator for the data
# In this way we make it continuous and do not need to carry u_d around.
from scipy import interpolate

u_d_fun=interpolate.interp1d(t_d,u_d,kind='cubic')

#%% Ingredients for the backward evaluation
# 1. simple forward pass + integration (function )
# 2. backward dynamical system (back_lambda) 
# 3. Backward integration; take in input the full history (Adjoint_Solver)
# 4. Compute J(p) and dJ_dp (J_Full_Back)

# General variables
n_t=1000; t0,tf=0,15; times=np.linspace(t0,tf,n_t)    


# This function implements the forward integration:
def forward(p):
    a,b,c,d=p[0],p[1],p[2],p[3]
    # We integrate with u0 and t: 
    u0=[1, 2]    
    # Here's the call to the ODE solver (RK4-5 by default)
    sol = solve_ivp(f,[t0,tf], u0, args=(a, b, c, d), t_eval=times,
                    dense_output=True,method='RK23')  
    # We provide the output as a callable function
    u=sol.sol(times)
    return u

# This function implements the ODE's function for lambda
def backward_lambda_AUG(t,lam_a,p):
    # unfold parameters... 
    p_1,p_2,p_3,p_4=p[0],p[1],p[2],p[3]  
    # ...the state variables
    u_1,u_2=lam_a[0],lam_a[1]
    # ...the adjoint variables
    lam=lam_a[2],lam_a[3]
    # ... and the available data
    u_1d,u_2d=u_d_fun(t)
    # unfold the state using the provided callable
    # Jacobian df_du:
    df_du=np.array([[p_1-p_2*u_2, -p_2*u_1],
                     [p_4*u_2, -p_3+p_4*u_1]])     
    # Gradient dL_du: we interpolate on the available data in time
    dL_du=2*np.array([u_1-u_1d,u_2-u_2d])
    
    # Advance the states
    u_1_dot=p_1*u_1 - p_2*u_1*u_2
    u_2_dot=-p_3*u_2 + p_4*u_1*u_2
    # Advance the adjoint states
    lam_dot =-df_du.T.dot(lam)-dL_du.T
    
    dlam_a_dt=np.hstack((u_1_dot,u_2_dot,lam_dot))
    
    return dlam_a_dt

# This function solves the adjoint problem for lambda
def Adjoint_Solver(p):
    # solve the forward problem for u(T)
    u=forward(p); u_T=u[:,-1]
    # Set the terminal conditions
    lam_a_T=np.hstack((u_T,0,0))
    # Here's the call to the ODE solver (RK4-5 by default)
    sol = solve_ivp(backward_lambda_AUG, 
                    [tf,t0], 
                    lam_a_T, args=(p,), t_eval=np.flip(times),
                    dense_output=True,method='RK23')  
    # Provide the output as a callable for lambda
    lam=sol.sol(np.flip(times))   
    return lam


#%% This implements the backward problem for J and dJ_dp
def J_Full_Back(p,u_d):
    # using the callabl u_fun, solve the backward pass    
    lam=Adjoint_Solver(p)
    # unwrap states, adjoint states and data along t=times
    u_1,u_2,lam_1,lam_2=lam; 
    u_1d,u_2d=u_d; u_1d=np.flip(u_1d);u_2d=np.flip(u_2d)
        
    # Compute the J(p) and dJ_
    integrand_J=np.zeros(n_t); integrand_dJ_dp=np.zeros((4,n_t))
    for k in range(n_t):
      integrand_J[k]= (u_1[k]-u_1d[k])**2+(u_2[k]-u_2d[k])**2 
      df_dp=np.array([[u_1[k], -u_1[k]*u_2[k], 0, 0],
                      [0,      0, -u_2[k],u_1[k]*u_2[k]]])
      integrand_dJ_dp[:,k]=np.array([lam_1[k],lam_2[k]]).T.dot(df_dp)
    # Finally compute J and dJ_dp
    J=np.sum(integrand_J)/n_t
    dJ_dp=np.sum(integrand_dJ_dp,axis=1)/n_t
    print('Cost F: {0:2f}'.format(J))
    return (J , dJ_dp)

(J_f , dJ_dp)=J_Full_Back(p_star,u_d)
print(' Evaluation with the backward approach on p='+str(p_star))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))

# This is hte result fo J and for p_0
(J_f , dJ_dp)=J_Full_Back(p0,u_d)
print(' Evaluation with the backward approach on p='+str(p0))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))
      

# This is hte result fo J and for p_2
p_2=np.array([1,2,3,3])
(J_f , dJ_dp)=J_Full_Back(p_2,u_d)

print(' Evaluation with the backward approach on p='+str(p_2))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))


#%% Optimization using Backward Gradient 
start_time = time.time()
solution = optimize.minimize(J_Full_Back,p0,method='BFGS',jac=True,args=(u_d,),
                              options={'disp': True, 'ftol':1e-19, 'maxiter':1000})
print("Finished in %s s" % (time.time() - start_time))   

p_final= solution.x # Optimal solution
print(p_final)

# Check out the results
u=forward(p_final); u_1,u_2=u


# Here's the plot of the solution    
fig, ax = plt.subplots(figsize=(6,4)) # Create Signal Noisy 
plt.plot(times,u_d[0,:])
plt.plot(times,u_d[1,:])
# The attempt from BFGS:
plt.plot(times,u_1,'k-',linewidth=2)
plt.plot(times,u_2,'k-',linewidth=2)

plt.xlabel('$t$',fontsize=12)
# plt.legend(shadow=True,fontsize=17)
#plt.title('Regression')
plt.tight_layout()
plt.show()
Name='BFGS_A_backward_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
# plt.close(fig)


fig, ax = plt.subplots(figsize=(4,4)) 
plt.plot(u_d[0,:],u_d[1,:],'r')
plt.plot(u_1,u_2,'k')

plt.xlabel('$u_1$',fontsize=14)
plt.ylabel('$u_2$',fontsize=14)
plt.tight_layout()
plt.show()
Name='BFGS_Backward_Orbit_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
plt.close(fig)







