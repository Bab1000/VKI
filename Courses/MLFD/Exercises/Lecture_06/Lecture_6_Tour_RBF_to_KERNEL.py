# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:51:05 2025

@author: mendez
"""

# This file provides a tour from ensamble learning in RBFs
# to Kernel Methods

# Import essential things
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

# plot customization (a matter of taste)

# Configuration for plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


#%% Section 0: Create the usual dataset
x1 = np.linspace(0, 4.3, 200, endpoint=True)
x2 = np.linspace(4.8, 10, 200, endpoint=True)
x_data=np.concatenate((x1,x2))
# Create the deterministic part
y_clean= 3*x_data+(x_data/100)**3+4*np.sin(3/2*np.pi*x_data)
# Add (a seeded) stochastic part
np.random.seed(0)
y=y_clean+1*np.random.randn(len(x_data))
# Introduce some outliers in x=2 and x=8
G1=10*np.exp(-(x_data-2)**2/0.005)*np.random.randn(len(x_data))
G2=15*np.exp(-(x_data-8)**2/0.005)*np.random.randn(len(x_data))
y_data=y+G1+G2

# Plot the data:

fig, ax = plt.subplots(figsize=(5, 3))
plt.scatter(x_data,y_data,c='white',
            marker='o',edgecolor='black',
            s=10,label='Data')
ax.set_xlabel('x data',fontsize=16)
ax.set_ylabel('y data ',fontsize=16)
plt.tight_layout()
plt.savefig('Sampled_data.png',dpi=200)


#%% Section 1: Define the functions for ensamble learning 

# This is the Gaussian RBF
def Gauss_RBF(x,x_r=0,c_r=0.1):
    d=x-x_r # Get distance
    phi_r=np.exp(-c_r**2*d**2)
    return phi_r

# Check the RBF
c_r=1 # Pick a shape factor
# Compute the RBF phi(x; x_c_r,c_r)
phi_0=Gauss_RBF(x_data,4,c_r=c_r)
# Here is how it looks like:
plt.plot(x_data,phi_0,'ko:')


# Prepare the Basis matrix using 100 equally spaced RBFs
n_b=100 # Define the number of bases
# Define where to put them:
x_b=np.linspace(x_data.min(),x_data.max(),n_b)

def PHI_Gauss_X(x_in, x_b, c_r=0.8):
 n_x=np.size(x_in)
 Phi_X=np.zeros((n_x,n_b+1)) # Initialize Basis Matrix on x
 # Add a constant and a linear term
 Phi_X[:,0]=x_in   # <-----------the linear term!
 # Loop to prepare the basis matrices (inefficient)
 for j in range(0,n_b):
  # Prepare all the terms in the basis
  Phi_X[:,j+1]=Gauss_RBF(x_in,x_r=x_b[j],c_r=c_r)
 return Phi_X

# Construct the basis function
Phi_X_s=PHI_Gauss_X(x_data, x_b, c_r=1)
# Plot them !
# plt.plot(x_data,Phi_X_s[:,0],'ko')

#% BOOTSTRAPPING Function
def Boot_strap_RBF_Train(x,y,x_b,c_r=0.6,alpha=0.001,n_e=500,tp=0.3):
  '''
  Bootstrap function that will traine n_e RBF models using Gaussian RBFs
  located at x_b, having shape factor c_r, using Ridge regression with alpha
  regularization and keeping a tp fraction for testing.

  '''
  J_i=np.zeros(n_e) # in sample error of the population
  J_o=np.zeros(n_e) # out of sample error of the population
  n_b=len(x_b)
  w_e=np.zeros((n_b+1,n_e)) # Distribution of weights

  # Loop over the ensamble
  for j in range(n_e):
         # Split the
     xs, xss, ys, yss = train_test_split(x,y, test_size=tp)
     # construct phi_x_s
     Phi_x_s=PHI_Gauss_X(xs, x_b, c_r=c_r)
     Phi_x_ss=PHI_Gauss_X(xss, x_b, c_r=c_r)

     # compute w
     H=Phi_x_s.T@Phi_x_s+alpha*np.identity(n_b+1)
     # Train Model
     w_s=np.linalg.inv(H).dot(Phi_x_s.T).dot(ys)

     # Assign vectors to the distributions
     w_e[:,j]=w_s
     # Make in-sample prediction---------------------------
     y_p_s=Phi_x_s.dot(w_s)
     # In-sample error
     J_i[j]= 1/(len(xs)-1)*np.linalg.norm(y_p_s-ys)**2
     # Make out-of sample prediction (and errors)
     y_p_ss=Phi_x_ss.dot(w_s)
     # Out of sample error
     J_o[j]=1/(len(xss)-1)*np.linalg.norm(y_p_ss-yss)**2
     # Fill the population matrix

  return J_i, J_o, w_e


def Ensamble_RBF(xg,n_b,x_b,c_r,w_e,sigma_y):
   '''
   Make enxamble prediction of models of RBF models using Gaussian RBFs
   located at x_b, having shape factor c_r from an ensamble of weights w_e
   on data that is estimated to have a random noise component sigma_y
   '''
   n_e=np.shape(w_e)[1] # get the n_e from the weight population
   n_p=len(xg) # number of points where predictions are required
   y_pop=np.zeros((n_p,n_e)) # Prepare the population of predictions

   # Prepare the Phi matrix on xg:
   Phi_X_ss=PHI_Gauss_X(xg, x_b, c_r=c_r)

   for e in range(n_e):
    y_pop[:,e]= Phi_X_ss.dot(w_e[:,e])

   # Get statistics over the ensamble----
   # Mean prediction
    y_e=np.mean(y_pop,axis=1)
    # the ensamble variance:
    Var_Y_model=np.std(y_pop,1)**2
    # So the global uncertainty, considering the aleatoric is:
    Unc_y=np.sqrt(Var_Y_model+sigma_y)
   return y_e, Unc_y

#-------------------------- PLAY WITH ME ------------------------------
# Let's see how to use these functions:
J_i, J_o, w_e=Boot_strap_RBF_Train(x_data,y_data,x_b,c_r=1,alpha=0.000001,n_e=100,tp=0.3)
sigma_y_estimate=np.sqrt(np.mean(J_i)) # Estimate the data variance

# Let's make prediction on 200 new points
x_ss=np.linspace(x_data.min(),x_data.max(),2000)
# Prediction of the ensamble
y_e, Unc_y=Ensamble_RBF(x_ss,n_b,x_b,c_r=1,w_e=w_e,sigma_y=sigma_y_estimate)


# Plot predictions
fig, ax = plt.subplots(figsize=(5, 3))
plt.scatter(x_data,y_data,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(x_ss,y_e,'r--',linewidth=2,label='$Ensamble Mean$')
plt.fill_between(x_ss, y_e + 1.96*Unc_y,
                 y_e - 1.96*Unc_y, alpha=0.3)

ax.set_xlabel('x',fontsize=12)
ax.set_ylabel('y',fontsize=12)
Name='Method_1_RBF_Ensamble_Predictions.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)



#%% Section 2: Ridge Regression + Weight Regression + Gaussian Propagation
# We fit a multivariate Gaussian on the population of weights
# Let us introduce the Gaussian Estimator (see slides ! )
def est_mu_Sigma(X):
   '''
   Estimate the mean and covariance from a population of vectors X (n_x,n_s)
   '''
   # the number of samples is
   n_s=X.shape[1]
   # The mean (as a vector) is
   mu_e=np.mean(X,axis=1).reshape(-1,1)
   # Mean shifted population
   X_s=X-np.repeat(mu_e,n_s,axis=1)
   Sigma_e=1/(n_s-1)*X_s.dot(X_s.T)
   return mu_e, Sigma_e

# This would be the estimator for the weight population:
mu_w,Sigma_w=est_mu_Sigma(w_e)

# Prediction on the training data to evaluate sigma_epsilon
Phi_X_s=PHI_Gauss_X(x_data, x_b, c_r=1)
y_data_s=Phi_X_s.dot(mu_w)[:,0]
sigma_y_s=np.std(y_data_s-y_data)

# The propagation on the unseen data woule then be:
Phi_X_ss=PHI_Gauss_X(x_ss, x_b, c_r=1)
# The linear propagation of Gaussians would give (see slides)
mu_y_w=Phi_X_ss.dot(mu_w)[:,0]; y_e=mu_y_w
Sigma_y_w=Phi_X_ss.dot(Sigma_w).dot(Phi_X_ss.T)

Sigma_Tot=Sigma_y_w+sigma_y_s**2*np.identity(len(mu_y_w))  # Check what happens if you remove the constant term (impact of unc!)
Unc_y = 1.96 * np.sqrt(np.diag(Sigma_Tot))


#%% FIGURE
fig, ax = plt.subplots(figsize=(5, 3))
plt.scatter(x_data,y_data,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(x_ss,y_e,'r--',linewidth=2,label='$Mean Posterior$')
plt.fill_between(x_ss, y_e + 1.96*Unc_y,y_e - 1.96*Unc_y, alpha=0.3)
ax.set_xlabel('x',fontsize=12)
ax.set_ylabel('y',fontsize=12)
Name='Method_2_Gauss_w_Prop.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)


#%% Section 3: Bayesian Ridge Regression

def Bayesian_RBF(x_s,y_s,x_ss,x_b,c_r=0.1,alpha=0.1):
   '''
   Bayesian RBF using training data (x_s,y_s), evaluated in x_ss, using a RBF
   Ridge regression with bases in x_b and shape factor c_r with regularization
   Ridge parameter alpha (see slides!)
   '''
   n_b=len(x_b) # number of RBFs elements

   #1. Compute the Basis Matrices
   PHI_X_s=PHI_Gauss_X(x_s, x_b, c_r=c_r)
   PHI_X_ss=PHI_Gauss_X(x_ss, x_b, c_r=c_r)
   #2. Prepare the regularized correlation matrix
   H=PHI_X_s.T@PHI_X_s+alpha*np.eye(n_b+1)
   Inv_H=np.linalg.inv(H)
   #3. Prepare the mu_w_y and Sigma_w_y from training
   mu_w_y=Inv_H.dot(PHI_X_s.T.dot(y_s)) # posterior on weights
   mu_y_s=PHI_X_s.dot(mu_w_y) # prediction on training data
   sigma_y=np.std(y_s.ravel()-mu_y_s) # estimate data variance on training
   # Posterior covariance
   Sigma_w_y=Inv_H*sigma_y**2

   # Prediction on x_ss
   mu_y=PHI_X_ss@mu_w_y
   Sigma_y=sigma_y**2*PHI_X_ss@Sigma_w_y@PHI_X_ss.T+\
             sigma_y**2*np.eye(len(x_ss))

   return mu_y, Sigma_y

#%% Bayesian Approach:
mu_y_w,Sigma_y_w=Bayesian_RBF(x_data,y_data,x_ss,x_b,c_r=1,alpha=0.001)
Unc = 1.96 * np.sqrt(np.diag(Sigma_y_w))


#%% FIGURE
fig, ax = plt.subplots(figsize=(5, 3))
plt.scatter(x_data,y_data,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(x_ss,mu_y_w,'r--',linewidth=2,label='$Mean Posterior$')
plt.fill_between(x_ss, mu_y_w + 1.96*Unc_y,mu_y_w - 1.96*Unc_y, alpha=0.3)
ax.set_xlabel('x',fontsize=12)
ax.set_ylabel('y',fontsize=12)
Name='Method_3_Bayesian_Ridge.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)

#%% Section 4: Kernelized Bayesian Ridge Regression

# We now introduce the Woodbury identity and Kernel functions!
from sklearn.metrics.pairwise import euclidean_distances

# Define Training/Testing
x_s, x_ss, y_s, y_ss = train_test_split(x_data,y_data, test_size=0.3)

# Completely new prediction: 500 new number
x_g=np.linspace(x_data.min(),x_data.max(),500)

def Kappa(X1,X2,l_c=1):
    sqdist = euclidean_distances(X1, X2)**2
    return  np.exp(-0.5 / l_c**2 * sqdist)

# We want to train only on x_s, then make predictions on x_g:
x_g=np.linspace(x_data.min(),x_data.max(),500)

# It is left to you to do the validation on (x_ss,y_ss)! 


def Kernel_Bayesian_RBF(x_s,y_s,x_ss,L_C,alpha=0.1):
   '''
   Kernelized Ridge Regression using training data (x_s,y_s), evaluated in x_ss, 
   using a RBF Kernel with length scale l_C and a regularization alpha (see slides!)
   '''
   # 1. Prepare the Kernel Matrices
   K_s_s=Kappa(x_s.reshape(-1,1),x_s.reshape(-1,1),l_c=L_C)
   K_ss_s=Kappa(x_ss.reshape(-1,1),x_s.reshape(-1,1),l_c=L_C)
   K_ss_ss=Kappa(x_ss.reshape(-1,1),x_ss.reshape(-1,1),l_c=L_C)

   # Prepare the inverted matrix (very inefficient!)
   H_u_s=K_s_s+alpha*np.eye(len(x_s)); Inv_H_u_s=np.linalg.inv(H_u_s)
   
   # 2. Compute the Mean predictions and estimate sigma_y 
   alpha_V=Inv_H_u_s.dot(y_s)
   mu_y_s=K_s_s.dot(alpha_V) # in sample prediction
   mu_y_ss=K_ss_s.dot(alpha_V) # out of sample prediction
   # Estimate sigma_y:
   sigma_y=np.std(mu_y_s-y_s)    
   
   # 3. Compute the Sigma_y:
   #efficient matrix multiplication
   block=np.linalg.multi_dot([K_ss_s,Inv_H_u_s,K_ss_s.T]) #
   Sigma_y_ss=sigma_y**2/alpha*K_ss_ss-sigma_y**2/alpha*block \
             + sigma_y**2*np.eye(len(x_ss))

   return mu_y_ss, Sigma_y_ss


# THE ONLY HYPER-PARAMETERS!
l_C=0.5 ; alpha=0.001

mu_y_ss, Sigma_y_ss=Kernel_Bayesian_RBF(x_s,y_s,x_g,l_C,alpha=0.1)

Unc_y = 1.96 * np.sqrt(np.diag(Sigma_y_ss))

# FIGURE
fig, ax = plt.subplots(figsize=(5, 3))
plt.scatter(x_data,y_data,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(x_g,mu_y_ss,'r--',linewidth=2,label='$Mean Posterior$')
plt.fill_between(x_g, mu_y_ss + 1.96*Unc_y,mu_y_ss - 1.96*Unc_y, alpha=0.3)
ax.set_xlabel('x',fontsize=12)
ax.set_ylabel('y',fontsize=12)
Name='Method_4_Kernel_RBF.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)















