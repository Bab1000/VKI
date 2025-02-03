# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:44:06 2025

@author: Miguel Alfonso Mendez
"""

# Part 1: We implement a linear model.

import numpy as np
import matplotlib.pyplot as plt
import pdb

# The following is an optional customization of the plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

# Step 1: Generate and plot the data

n_p=100
np.random.seed(10)
x_s=np.random.uniform(0,10,n_p)
y_s=2*x_s+2+np.random.normal(loc=0,scale=10,size=len(x_s))+x_s**2

# Show the results
fig, ax = plt.subplots(figsize=(5, 3)) # This creates the figure
plt.scatter(x_s,y_s,c='black',marker='o',edgecolor='black',s=16)
ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Exercise_1_data.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 

# We save the data for later use:
np.savez('Data_Ex_1',y_s=y_s,x_s=x_s)

    
    
# Create functions for regression

def line_s_fit(x,y):
    # take input vectors x,y and fit a line
    # with coefficients w_0 and w_1 
    n_p=len(x); X=np.zeros((n_p,2)) # get number of points and initialize X

    X[:,0] = 1
    X[:,0] = x

    H = X.T@(X)
    w = np.linalg.inv(H).dot(X.T).dot(y)
    return w, H    

def Cubic_s_fit(x,y):
    # take input vectors x,y and fit a quadratic
    # with coefficients w_0, w_1, w_2 and w_3
    n_p=len(x); X=np.zeros((n_p,4))

    X[:,0] = 1
    X[:,0] = x
    X[:,0] = x**2
    X[:,0] = x**3

    H = X.T@(X)
    w = np.linalg.inv(H).dot(X.T).dot(y)

    return w, H    



# Create a function for the general model of order n_O
def General_s_fit(x,y,n_O):
   # take input vectors x,y and fit a polynomial
   # with coefficients w_0, w_1,...w_{n_O}
    n_p=len(x)
    X=np.zeros((n_p,n_O+1))
    for j in range(n_O+1):
     X[:,j]=x**j   
    # Compute the Hessian
    H=X.T@(X)
    w=np.linalg.inv(H).dot(X.T).dot(y)
    return w,H



# Variant with the Cholesky decomposition:
    
# Here is an example
from scipy.linalg import cholesky, cho_solve
# Define the symmetric positive definite matrix A 
# and the right-hand side vector b
A = np.array([[4, 2, 2], 
              [2, 6, 2], 
              [2, 2, 5]])
b = np.array([12, 18, 17])

# Perform Cholesky decomposition
L = cholesky(A, lower=True)  # `lower=True` computes the lower triangular L

# Solve the system using the Cholesky factorization
# Forward substitution for Ly = b, followed by backward substitution for L^T x = y
x = cho_solve((L, True), b)  # (L, True) indicates L is lower triangular

print("Solution x:", x)

def General_s_fit_Chol(x,y,n_O):
   # take input vectors x,y and fit a polynomial
   # with coefficients w_0, w_1,...w_{n_O}, assuming that
   # the data points have a variance sigma_y
    n_p=len(x)
    X=np.zeros((n_p,n_O+1))
    for j in range(n_O+1):
     X[:,j]=x**j   
    # Compute the Hessian and the RHS
    H=X.T@(X); b=X.T.dot(y)
    L= cholesky(H, lower=True)
    w= cho_solve((L,True),b)
    return w, H


# Variant with distributed uncertainties:

def General_s_fit_Chol_w(x,y,sigma_y,n_O):
   # take input vectors x,y and fit a polynomial
   # with coefficients w_0, w_1,...w_{n_O}, assuming that
   # the data points have a variance sigma_y
    n_p=len(x)
    X=np.zeros((n_p,n_O+1))
    S=np.diag(1/sigma_y**2)
    for j in range(n_O+1):
     X[:,j]=x**j   
    # Compute the Hessian and the RHS
    H=X.T@S@(X); b = X.T@S.dot(y)
    L=cholesky(H,lower=True) # Cholesky Deco
    w=cho_solve((L,True),b)
    return w, H



# Compare the results with numpy's polyfit

# 1. Linear case. Your implementation:
print(' \nValidation for the linear model\n')
w_lin,H_lin=General_s_fit(x_s,y_s,1)
print('w0, w1 from your implementation:{:.3f}, {:.3f}'.format(w_lin[0],w_lin[1]))
# linear, Python's implementation
w_lin_np=np.polyfit(x_s,y_s,1)
print('w0, w1 from numpy:{:.3f}, {:.3f}'.format(w_lin_np[1],w_lin_np[0]))

# 2. Quadratic case. Your implementation:
print(' \nValidation for the Quadratic model\n')
w_quad,H_quad=General_s_fit(x_s,y_s,2)
print('w0, w1 ,w2 from your implementation:{:.3f}, {:.3f}, {:.3f}'
      .format(w_quad[0],w_quad[1],w_quad[2]))
# quadratic, Python's implementation
w_quad_np=np.polyfit(x_s,y_s,2)
print('w0, w1, w2 from numpy:{:.3f}, {:.3f}, {:.3f}'.
      format(w_quad_np[2],w_quad_np[1],w_quad_np[0]))


# 3. Cubic case. Your implementation:
print(' \nValidation for the cubic model\n')
w_cub,H_cub=General_s_fit(x_s, y_s,3)
print('w0, w1 ,w2, w3 from your implementation:{:.3f}, {:.3f}, {:.3f}, {:.3f}'
      .format(w_cub[0],w_cub[1],w_cub[2],w_cub[3]))
# quadratic, Python's implementation
w_cub_np=np.polyfit(x_s,y_s,3)
print('w0, w1 ,w2, w3 from numpy :{:.3f}, {:.3f}, {:.3f}, {:.3f}'
      .format(w_cub_np[3],w_cub_np[2],w_cub_np[1],w_cub_np[0]))


# Final check to close this step: What is the cond number for these three?

print('\nCondition number for H in model 1:{:.4e}'.format(np.linalg.cond(H_lin)))

print('\nCondition number for H in model 2:{:.4e}'.format(np.linalg.cond(H_quad)))

print('\nCondition number for H in model 3:{:.4e}'.format(np.linalg.cond(H_cub)))


# Make prediction on the available data
def Poly_model_Pred(x,w):
    # Make predictions of x using the polynomial model with weights w
    # Check the size of our vectors 
    n_p=len(x); n_b=len(w)
    # Initialize X:
    X=np.zeros((n_p,n_b))
    # Build X    
    for j in range(n_b):
     X[:,j]=x**j    
    y=X.dot(w)
    return y

y_lin=Poly_model_Pred(x_s,w_lin)
y_quad=Poly_model_Pred(x_s,w_quad)
y_cub=Poly_model_Pred(x_s,w_cub)


# Show the results
fig, ax = plt.subplots(figsize=(5, 3)) # This creates the figure
plt.scatter(x_s,y_s,c='white',marker='o',edgecolor='black',s=16)
plt.plot(x_s,y_lin,'rs',label='lin')
plt.plot(x_s,y_quad,'bo',label='quad')
plt.plot(x_s,y_cub,'kv',label='cub')
plt.legend(fontsize=14)
ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Exercise_1_data_with_fits.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 


# Predictions with weighted least squares

sigma_y=1/2*x_s**2 + 5

# linear model with constant uncertainties
w_lin,H_lin=General_s_fit(x_s,y_s,1)
y_lin=Poly_model_Pred(x_s,w_lin)

# linear model with varying uncertainties
w_lin_W,H_lin_W=General_s_fit_Chol_w(x_s,y_s,sigma_y,1)
y_lin_W=Poly_model_Pred(x_s,w_lin_W)


# Show the results
fig, ax = plt.subplots(figsize=(5, 3)) # This creates the figure
plt.errorbar(x_s, y_s, yerr=1.96*sigma_y, xerr=0, fmt='o', capsize=5, label='Data with error bars')
plt.plot(x_s,y_lin,'rs',label='lin')
plt.plot(x_s,y_lin_W,'bo',label='lin w')
plt.legend(fontsize=14)
ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Exercise_1_data_with_fits_varying_sigma.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 


