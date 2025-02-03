# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 20:28:40 2023

@author: mendez
"""

import numpy as np

import matplotlib.pyplot as plt

# The following is an optional customization of the plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


# We save the data for later use:
data = np.load('Data_Ex_1.npz')
y_s=data['y_s']; x_s=data['x_s']
n_p=len(x_s)
    

# Create a cost function

def cost(w,X,y):
    n_p = np.shape(y)[0]
    J = 1/n_p*(np.linalg.norm(y-X.dot(w)))**2
  
    return J

# Create a gradient function

def grad(w,X,y):
    n_p = np.shape(y)[0]
    Nabla_J = 2/n_p*X.T.dot(X.dot(w) - y)

    return Nabla_J



def Batch_GD(cost,grad,w0,x,y,eta_0,decay,n_epochs,n_batch):
    # Implementation of the SGD training for a polynomial model with
    # cost function 'cost', having gradieng 'grad' and 
    # starting from initial guess of weights 'w0'.
    # The training data are (x,y). The learning follows a decay law, governed 
    #  by the parameter 'decay', and starting from learning rate eta_0.
    # The training is carried out for n_epochs epochs and using
    # batches of size n_batch.
    
    # get the number of features and the number of points
    n_b=len(w0); n_p=len(x)
    # Prepare the complete feature matrix
    X=np.zeros((n_p,n_b))
    # Build X    
    for j in range(n_b):
     X[:,j]=x**j
       
    # Prepare the loop per epoch
    n_iter=n_epochs*n_p//n_batch     
    # Current estimate of w
    w=w0
    # Initialize the weight evolution and the error evolution
    Err_SGD=np.zeros(n_iter); #Err_SGD[0]=cost(w,X,y)
    w_evolution=np.zeros((n_b,n_iter)); #w_evolution[:,0]=w0
    
    for j in range(n_iter):      
      # Select randomly some data points for the batch
      # Note that replace=False means that there is no repetition
      Indices=np.random.choice(n_p, n_batch,replace=False)
      # Construct the matrix X_b
      X_b=X[Indices,:]; y_b=y[Indices]  
      # Get the current cost
      Err_SGD[j]=cost(w,X_b,y_b)
      #Get the gradient
      Nabla_J_w=grad(w,X_b,y_b)
      # Compute the learning rate
      eta=eta_0/(1+decay*j)
      # Weght update
      w=w-eta*Nabla_J_w
      # Store the result in the history
      w_evolution[:,j]=w
      # Message
      Mex='Iteration: {:d}; Epoch: {:d}; Cost: {:.3f}; Grad_abs: {:.3f}'\
          .format(j,(j*n_batch)//n_p,Err_SGD[j],np.linalg.norm(Nabla_J_w))
      print(Mex)    
      
    w_opt=w # Final result on the weight  
    return w_opt, w_evolution, Err_SGD




#%% Test on linear model with varioys batch sizes

# Take first order
w0=np.array([0, 0 ])


w_opt, w_evolution_1, Err_SGD_1=Batch_GD(cost,grad,w0,x_s,y_s,eta_0=10e-4,decay=0,
                                     n_epochs=10000,n_batch=100)


w_opt, w_evolution_2, Err_SGD_2=Batch_GD(cost,grad,w0,x_s,y_s,eta_0=10e-4,decay=0,
                                     n_epochs=10000//20,n_batch=5)

# We define the X here to call 'cost' on the grid of the parameter space
n_b=len(w0); n_p=len(x_s)
# Prepare the complete feature matrix
X=np.zeros((n_p,n_b))
# Build X    
for j in range(n_b):
   X[:,j]=x_s**j



# Plot the result in the weight space:
w_1_g=np.linspace(-30,80,200)
w_2_g=np.linspace(-10,20,200)
    
w1g,w2g=np.meshgrid(w_1_g,w_2_g)
J=np.zeros_like((w1g))    
# Evaluate J in all of these points!
for i in range(len(w_1_g)):
  for j in range(len(w_2_g)):
    J[i,j]=cost(np.array([w1g[i,j],w2g[i,j]]),X,y_s)         

fig, ax = plt.subplots(figsize=(4, 4)) # This creates the figure
plt.contourf(w1g,w2g,J,alpha=0.7)
plt.plot(w_evolution_1[0,:],w_evolution_1[1,:],'r-',linewidth=2,label='GD')
plt.plot(w_evolution_2[0,:],w_evolution_2[1,:],'b--',linewidth=2,label='SGD')
ax.set_xlabel('$w_0$',fontsize=14)
ax.set_ylabel('$w_1$',fontsize=14)
plt.legend()
plt.tight_layout()
Name='Trajectory.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)
plt.close("all")
    
    
fig, ax = plt.subplots(figsize=(6, 3)) # This creates the figure
plt.plot(Err_SGD_2,'b--',linewidth=2,label='SGD')
plt.plot(Err_SGD_1,'r-',linewidth=2,label='GD')

ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlabel('Iterations',fontsize=14)
ax.set_ylabel('$J(\mathbf{w})$',fontsize=14)
plt.legend()
plt.tight_layout()
Name='convergence_GD_vs_SGD.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)
plt.close("all")
    
    













