# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:55:56 2023

@author: mendez
"""

# We perform the Assimilation as before, but using a forward approach to 
# the gradient evaluation

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
t=data['t']
u_d=data['u_d']

#%% Ingredients for the forward evaluation
# 1. Dynamical system that has both both the states AND the sensitivities (function full_S )
# 2. Forward pass that advances both (forward_full)
# 3. Compute J(p) and dJ_dp


#%% Full system with forward evaluation of sensitivities
# This function implements the augmented system (with sensitivities) for the
# Lotka- Volterra system 
def full_S(t,u_a,p):
    # unfold the set of parameters:
    p_1,p_2,p_3,p_4=p[0],p[1],p[2],p[3]    
    # unfold variables and sensitivities (each of which is a vector of 2)
    u_1=u_a[0];  u_2=u_a[1]
    s_1=u_a[2:4] # sensitivity vector with respect to p1
    s_2=u_a[4:6] # sensitivity vector with respect to p2
    s_3=u_a[6:8] # sensitivity vector with respect to p3
    s_4=u_a[8:10] # sensitivity vector with respect to p4
    
    # Derivative of f with respect to u
    df_du=np.array([[p_1-p_2*u_2, -p_2*u_1],
                    [p_4*u_2, -p_3+p_4*u_1]])    
    # Derivative of f with respect to p
    df_dp=np.array([[u_1, -u_1*u_2, 0, 0],
                    [0, 0,-u_2,u_1*u_2]])
    
    # Advance the system 
    u_1_dot=p_1*u_1 - p_2*u_1*u_2; u_2_dot=-p_3*u_2 + p_4*u_1*u_2
    
    # advance the sensitivities
    s_1_dot=df_du.dot(s_1)+df_dp[:,0]; s_2_dot=df_du.dot(s_2)+df_dp[:,1]
    s_3_dot=df_du.dot(s_3)+df_dp[:,2]; s_4_dot=df_du.dot(s_4)+df_dp[:,3]
    
    # Pile up the whole vector:
    u_a_dot=np.hstack((u_1_dot,u_2_dot,s_1_dot,s_2_dot,s_3_dot,s_4_dot))
    
    return u_a_dot

#%% Forward, full system ODE integration 
# Here we have the full forward evaluation from the parameters p

def forward_full(p):
    # We integrate with u0 and t: 
    u_a_0=[1, 2,0,0,0,0,0,0,0,0]; t0,tf=0,15; t=np.linspace(t0,tf,1000)    
    # Here's the call to the ODE solver (RK4-5 by default)
    sol = solve_ivp(full_S, 
                    [t0,tf], 
                    u_a_0, args=(p,), t_eval=t,
                    dense_output=True,method='RK23')  
    # Here's the solution of the system:
    u_a=sol.sol(t)  
    return u_a     

#%% This is the cost function that gives both J and dJ_dp
def J_FULL_F(p,u_d):
    # Do a forward pass with the parameters p
    u_a=forward_full(p); n_s=np.shape(u_a)[1]
    # Re-assign the variables...
    u_1=u_a[0,:]; u_2=u_a[1,:]
    # ... and the sensitivities
    du1_dp1=u_a[2,:]; du2_dp1=u_a[3,:]
    du1_dp2=u_a[4,:]; du2_dp2=u_a[5,:]
    du1_dp3=u_a[6,:]; du2_dp3=u_a[7,:]
    du1_dp4=u_a[8,:]; du2_dp4=u_a[9,:]
    # Compute the cost function:
    L=(u_1-u_d[0,:])**2+(u_2-u_d[1,:])**2; J=np.sum(L)/n_s    
    # ... and compute the gradient. We proceed term by term:
    dJ_dp1=0; dJ_dp2=0;dJ_dp3=0;dJ_dp4=0
    #.... which we compute via summation
    for j in range(len(u_1)):
      delta_1_j=(u_1[j]-u_d[0,j]); delta_2_j=(u_2[j]-u_d[1,j])
      
      dJ_dp1+=2*delta_1_j*du1_dp1[j]+2*delta_2_j*du2_dp1[j]
      dJ_dp2+=2*delta_1_j*du1_dp2[j]+2*delta_2_j*du2_dp2[j]
      dJ_dp3+=2*delta_1_j*du1_dp3[j]+2*delta_2_j*du2_dp3[j]
      dJ_dp4+=2*delta_1_j*du1_dp4[j]+2*delta_2_j*du2_dp4[j]
      
    dJ_dp=np.array([dJ_dp1,dJ_dp2,dJ_dp3,dJ_dp4])/n_s   
    print('Cost F: {0:2f}'.format(J))
    return (J , dJ_dp)
    

(J_f , dJ_dp)=J_FULL_F(p_star,u_d)
print(' Evaluation with the Forward approach on p='+str(p_star))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))

# This is hte result fo J and for p_0
(J_f , dJ_dp)=J_FULL_F(p0,u_d)

print(' Evaluation with the Forward approach on p='+str(p0))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))


# This is the result fo J and for p_2
p_2=np.array([1,2,3,3])
(J_f , dJ_dp)=J_FULL_F(p_2,u_d)

print(' Evaluation with the Forward approach on p='+str(p_2))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))


      
#%% Optimization using Forward Sensitivities
start_time = time.time()
solution = optimize.minimize(J_FULL_F,p0,method='BFGS',jac=True,args=(u_d,),
                             options={'disp': True})
print("Finished in %s s" % (time.time() - start_time))   
p_final= solution.x # Optimal solution
print(p_final)

# Check the results:
u_a_BFGS=forward_full(p_final)

t0,tf=0,15; t=np.linspace(t0,tf,1000)  


# Here's the plot of the solution    
fig, ax = plt.subplots(figsize=(6,4)) # Create Signal Noisy 
plt.plot(t,u_d[0,:],label='$\mathbf{u}_1$')
plt.plot(t,u_d[1,:],label='$\mathbf{u}_2$')
# The attempt from BFGS:
plt.plot(t,u_a_BFGS[0,:],'k-',linewidth=2,label='$\mathbf{u}_1$ BFGS')
plt.plot(t,u_a_BFGS[1,:],'k-',linewidth=2,label='$\mathbf{u}_2$ BFGS')

plt.xlabel('$t$',fontsize=12)
# plt.legend(shadow=True,fontsize=17)
#plt.title('Regression')
plt.tight_layout()
plt.show()
Name='BFGS_A_Forward_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
# plt.close(fig)

fig, ax = plt.subplots(figsize=(4,4)) 
plt.plot(u_d[0,:],u_d[1,:],'r')
plt.plot(u_a_BFGS[0,:],u_a_BFGS[1,:],'k')

plt.xlabel('$u_1$',fontsize=14)
plt.ylabel('$u_2$',fontsize=14)
plt.tight_layout()
plt.show()
Name='BFGS_Forward_Orbit_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
plt.close(fig)





