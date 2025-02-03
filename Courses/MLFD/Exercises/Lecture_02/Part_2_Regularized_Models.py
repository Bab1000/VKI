# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 18:07:40 2023

@author: Miguel Alfonso Mendez
"""

# Part 1: We implement a linear model.

import numpy as np
import matplotlib.pyplot as plt


# Import the function for the prediction from the previous file
from Part_1_Linear_Model import Poly_model_Pred


# The following is an optional customization of the plots
plt.rc('text', usetex=False)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

# Step 1: Load the data

# We save the data for later use:
data = np.load('Data_Ex_1.npz')
y_s=data['y_s']; x_s=data['x_s']
n_p=len(x_s)
    
    
# Create functions for regularized regression


def General_s_fit_reg(x,y,n_O,k_l=1000):
   # take input vectors x,y and fit a polynomial
   # with coefficients w_0, w_1,...w_{n_O} and imposing
   # an upper limit k_l to the condition number
    n_p=len(x)
    X=np.zeros((n_p,n_O+1))
    for j in range(n_O+1):
     X[:,j]=x**j   
    # Compute the Hessian
    H=X.T@(X)
    # compute the eigenvalues, the largest (l_M) and the smallest (l_m)
    Lambd=np.linalg.eig(H); l_M=np.max(Lambd[0]); l_m=np.min(Lambd[0])
    # compute the alpha
    if l_M/l_m > k_l:
     alpha = (l_M - k_l*l_m)/(k_l - 1)
     print("Regularization applied !")

    else:
     alpha= 0
     
    # regularize the Hessian ( if needed!)
    H_p=H+alpha*np.identity(np.shape(H)[0])
    # solve for w:
    w=np.linalg.inv(H_p).dot(X.T).dot(y)
    return w, H_p, alpha    


# Compare the results with numpy's polyfit

K_Lim=2e5 # Define the regularization for all problems.

# 1. Linear case. Your implementation:
print(' \nRegularization of the linear model\n')
w_lin,H_p_lin, alpha_lin=General_s_fit_reg(x_s,y_s,1,k_l=K_Lim)
print('w0, w1 with k_l={:.1e} :{:.3f}, {:.3f}'.format(K_Lim,w_lin[0],w_lin[1]))
print('\n computed with alpha={:.3f}'.format(alpha_lin))
print('\nCondition number for H_p in model 1:{:.4e}'.format(np.linalg.cond(H_p_lin)))

print('------------------------------------------')
print('------------------------------------------')
print('------------------------------------------')

# 2. Cubic case. 
print(' \nRegularization of the cubic model\n')
w_cub,H_p_cub, alpha_cub=General_s_fit_reg(x_s,y_s,3,k_l=K_Lim)
print('w0, w1, w2, w3 with k_l={:.1e} :{:.3f}, {:.3f}, {:.3f}, {:.3f}'
      .format(K_Lim,w_cub[0],w_cub[1],w_cub[2],w_cub[3]))
print('\n computed with alpha={:.3f}'.format(alpha_cub))
print('\nCondition number for H_p in model 1:{:.4e}'.format(np.linalg.cond(H_p_cub)))




# Make prediction on the available data
y_lin=Poly_model_Pred(x_s,w_lin)
y_cub=Poly_model_Pred(x_s,w_cub)



# Show the results
fig, ax = plt.subplots(figsize=(5, 3)) # This creates the figure
plt.scatter(x_s,y_s,c='white',marker='o',edgecolor='black',s=16)
plt.plot(x_s,y_lin,'rs',label='lin')
plt.plot(x_s,y_cub,'kv',label='cub')
plt.legend(fontsize=14)

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Exercise_1_data_with_fits.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 














