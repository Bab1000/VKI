# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:57:06 2023

@author: mendez
"""

# In this exercise we seek to perform data assimilation on the Lotka-Volterra
# model. We here perform first a crude evaluation of the cost function gradient
# using finite differences and run a BFGS optimization with it

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import optimize
import time


# Configuration for plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
 

# Here is the definition of the system with the default parameters
a,b,c,d=1,1,3,1;

# Define the 'Good parameters':
p_star=np.array([a,b,c,d])


# Define the 'Initial Guesses' for the inverse problem
p0=np.array([1.2,1.4,3.4,1.2])

def f(t, u, a=1, b=1, c=3, d=1):
    u1, u2 = u
    return [a*u1 - b*u1*u2, -c*u2 + d*u1*u2]


#%% Generate the (noisy data)
# We can integrate this forward in time for t in [0,15], starting from 
u0=[1, 2]    
# and on a vector of time of 1000 points:
t0,tf=0,15    
t=np.linspace(t0,tf,1000)    
# Here's the call to the ODE solver ( RK4-5 by default)
sol = solve_ivp(f, 
                [t0,tf], 
                u0, args=(a,b,c,d), t_eval=t,
                dense_output=True)  
# Here's the solution of the system:
u=sol.sol(t)    

#%%%%%%%%%%%%%%%%%% Data Generation/ Plot ####################
u_d=u+np.random.rand(2,1000)*0.5
########################################################


# Here's the plot of the solution    
fig, ax = plt.subplots(figsize=(6,4)) # Create Signal Noisy 
u_N=u+np.random.rand(2,1000)*0.5
plt.plot(t,u_d[0,:],label='$\mathbf{u}_1$')
plt.plot(t,u_d[1,:],label='$\mathbf{u}_2$')
plt.xlabel('$t$',fontsize=14)
plt.legend(shadow=True)
plt.title('Solution of the Lotka-Volterra System')
plt.tight_layout()
plt.show()
Name='Solution_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
plt.close(fig)
    

fig, ax = plt.subplots(figsize=(4,4)) # Create Signal Noisy 
plt.plot(u_d[0,:],u_d[1,:])
plt.xlabel('$u_1$',fontsize=14)
plt.ylabel('$u_2$',fontsize=14)
plt.tight_layout()
plt.show()
Name='Orbit_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
plt.close(fig)


#%% Save the data for later use:

np.savez('Data_Ex_Assimilation',t=t,u_d=u_d)


#%% Definition of the cost function and its gradient

# This function implements the forward integration:
def forward(p):
    a,b,c,d=p[0],p[1],p[2],p[3]
    # We integrate with u0 and t: 
    u0=[1, 2]; t0,tf=0,15; t=np.linspace(t0,tf,1000)    
    # Here's the call to the ODE solver (RK4-5 by default)
    sol = solve_ivp(f, 
                    [t0,tf], 
                    u0, args=(a, b, c, d), t_eval=t,
                    dense_output=True,method='RK23')  
    # Here's the solution of the system:
    u=sol.sol(t)  
    return u 

# This is the problem's cost function 
def J(p,u_d):
    # Do a forward pass with the parameters p
    u=forward(p); n_s=np.shape(u)[1]    
    # Now define the cost function:
    L=(u[0,:]-u_d[0,:])**2+(u[1,:]-u_d[1,:])**2
    J=np.sum(L)/n_s 
    return J

# We can now implement the full cost function with FD gradient:

def J_FULL_Fin_D(p,u_d):    
    # evaluate the cost function    
    J_value=J(p,u_d) 
    # Evaluate the gradient using FD
    dJ_dp=optimize.approx_fprime(p,J,[1.5e-9,1.5e-9,1.5e-9,1.5e-9],u_d)  
    print('Cost F: {0:2f}'.format(J_value)) 
    return (J_value, dJ_dp)


# This is hte result fo J and for p_0
(J_f , dJ_dp)=J_FULL_Fin_D(p_star,u_d)

print(' Evaluation with the Fin D approach on p='+str(p_star))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))


# This is the result fo J and for p_0
(J_f , dJ_dp)=J_FULL_Fin_D(p0,u_d)

print(' Evaluation with the Fin D approach on p='+str(p0))
print('J='+ str(J_f))
print('dJ_dp='+str(dJ_dp))

    
#%% Here's a very crude optimization using numpy's BFGS with finite differences
start_time = time.time()
solution = minimize(J_FULL_Fin_D,p0,method='BFGS',args=(u_d,),
                    jac=True,
                    options={'disp': True})
print("Finished in %s s" % (time.time() - start_time))   


p_final= solution.x # Optimal solution
u_BFGS=forward(p_final)

# Here's the plot of the solution    
fig, ax = plt.subplots(figsize=(6,4)) 
plt.plot(t,u_N[0,:],label='$\mathbf{u}_1$')
plt.plot(t,u_N[1,:],label='$\mathbf{u}_2$')
# The attempt from BFGS:
plt.plot(t,u_BFGS[0,:],'k-',linewidth=2,label='$\mathbf{u}_1$ BFGS')
plt.plot(t,u_BFGS[1,:],'k-',linewidth=2,label='$\mathbf{u}_2$ BFGS')

plt.xlabel('$t$',fontsize=12)
# plt.legend(shadow=True,fontsize=17)
#plt.title('Regression')
plt.tight_layout()
plt.show()
Name='BFGS_A_FD_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
# plt.close(fig)



fig, ax = plt.subplots(figsize=(4,4)) 
plt.plot(u_d[0,:],u_d[1,:],'r')
plt.plot(u_BFGS[0,:],u_BFGS[1,:],'k')

plt.xlabel('$u_1$',fontsize=14)
plt.ylabel('$u_2$',fontsize=14)
plt.tight_layout()
plt.show()
Name='BFGS_Orbit_'+str(a)+'_'+str(b)+'_'+str(c)+'_'+str(d)+'.png'
plt.savefig(Name, dpi=300)      
plt.close(fig)
