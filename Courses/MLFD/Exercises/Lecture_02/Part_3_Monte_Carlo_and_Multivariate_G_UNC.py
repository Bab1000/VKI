# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 18:07:40 2023

@author: Miguel Alfonso Mendez
"""

# Part 1: We implement a linear model.

import numpy as np
import matplotlib.pyplot as plt

# The function we created before to train the quadratic model
from Part_2_Regularized_Models import General_s_fit_reg

# The function we created before to make prediction using a polynomial model
from Part_1_Linear_Model import Poly_model_Pred


# Useful function for the splitting into training and testing
from sklearn.model_selection import train_test_split


# The following is an optional customization of the plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

#%% Step 1: Train the ensemble model using boot-strapping

# We save the data for later use:
data = np.load('Data_Ex_1.npz')
y_s=data['y_s']; x_s=data['x_s']
n_p=len(x_s)
    

def Boot_strap_Poly_Train(x,y,n_O,k_l=1e3,n_e=500,tp=0.3):
    """
    Probabilistic polynomial model of order n_O, trained on the data x,y,
    with upper condition number k_l, using a population of n_e realization
    and a test-splitting with ratio tp
    """

    J_i=np.zeros(n_e) # in sample error of the population
    J_o=np.zeros(n_e) # out of sample error of the population
    w_e=np.zeros((n_O+1,n_e)) # Distribution of weights

   # Loop over the ensamble    
    for j in range(n_e):    
     # Split the    
     xs, xss, ys, yss = train_test_split(x,y, test_size=tp)    # s= training, ss = testing
     # Fit the weights on the training data (xs,ys)
     w_s, _, _  =General_s_fit_reg(xs,ys,n_O,k_l=k_l)
     # Assign vectors to the distributions
     w_e[:,j]=w_s 
     # Make in-sample prediction---------------------------
     y_p_s=Poly_model_Pred(xs,w_s); 
     # In-sample error
     J_i[j] = 1/len(xs)*np.linalg.norm(y_p_s-ys)**2
     # Make out-of sample prediction (and errors)
     y_p_ss=Poly_model_Pred(xss,w_s)
     # Out of sample error
     J_o[j] = 1/len(xss)*np.linalg.norm(y_p_ss-yss)**2
     # Fill the population matrix
         
    return J_i, J_o, w_e


# Train a probabilistic linear model on the available data
J_i_1, J_o_1, w_e_1=Boot_strap_Poly_Train(x_s,y_s,1,k_l=1e6,n_e=500,tp=0.3)


# Train a probabilistic cubic model on the available data
J_i_2, J_o_2, w_e_2=Boot_strap_Poly_Train(x_s,y_s,2,k_l=1e6,n_e=500,tp=0.3)

# Train a probabilistic cubic model on the available data
J_i_3, J_o_3, w_e_3=Boot_strap_Poly_Train(x_s,y_s,3,k_l=1e6,n_e=500,tp=0.3)


#% Plot the results:

# In-sample and out of sample errors for the linear model 
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.hist(J_o_1,100,label='$J_o$',alpha=0.8)
plt.hist(J_i_1,100,label='$J_i$',alpha=0.8)
ax.set_xlabel('$J_i$, $J_o$,',fontsize=14)
ax.set_ylabel('$p(J_i)$, $p(J_o)$',fontsize=14)
plt.legend()
plt.tight_layout()
Name='J_i_J_o_1.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)
plt.close("all")


# In-sample and out of sample errors for the linear model 
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.hist(J_o_2,100,label='$J_o$',alpha=0.8)
plt.hist(J_i_2,100,label='$J_i$',alpha=0.8)
ax.set_xlabel('$J_i$, $J_o$,',fontsize=14)
ax.set_ylabel('$p(J_i)$, $p(J_o)$',fontsize=14)
plt.legend()
plt.tight_layout()
Name='J_i_J_o_2.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)
plt.close("all")

# In-sample and out of sample errors for the cubic model 
fig, ax = plt.subplots(figsize=(5, 3))
plt.hist(J_o_3,100,label='$J_o$',alpha=0.8)
plt.hist(J_i_3,100,label='$J_i$',alpha=0.8)
ax.set_xlabel('$J_i$, $J_o$,',fontsize=14)
ax.set_ylabel('$p(J_i)$, $p(J_o)$',fontsize=14)
plt.legend()
plt.tight_layout()
Name='J_i_J_o_3.png'
plt.tight_layout()
plt.savefig(Name, dpi=300)
plt.close("all")


# Uncertainty evaluation on new points

xg=np.linspace(x_s.min(),x_s.max(),200)


def Ensamble_Pred(xg,w_e,J_i_mean):
    # Make prediction on the points xg for a polynomial model with 
    # ensamble weights w_e and mean in-sample error of J_i_mean
    n_Op1,n_e=np.shape(w_e); n_p=len(xg)
    # prepare the population of predictions in y:
    y_pop=np.zeros((n_p,n_e))    
    for j in range(n_e):   # loop over the enxamble  
     y_pop[:,j]= Poly_model_Pred(xg,w_e[:,j])

    # The mean prediction will be:
    y_e=np.mean(y_pop,axis=1)   
    # the ensamble std:
    Var_Y_model=np.std(y_pop,1)**2    
    # So the global uncertainty, considering the aleatoric is:
    Unc_y=np.sqrt(Var_Y_model+J_i_mean)
    
    return y_e, Unc_y
                               
                               
# Predictions of the ensamble of linear model:                              
y_e_1, Unc_y_1=Ensamble_Pred(xg,w_e_1,np.mean(J_i_1))

# Predictions of the ensamble of quadratic models:                              
y_e_2, Unc_y_2=Ensamble_Pred(xg,w_e_2,np.mean(J_i_2))

# Predictions of the ensamble of cubic models:                              
y_e_3, Unc_y_3=Ensamble_Pred(xg,w_e_3,np.mean(J_i_3))



# Plot linear model with Uncertainties
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.scatter(x_s,y_s,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(xg,y_e_1,'r--',linewidth=2)
plt.fill_between(xg, y_e_1 + 1.96*Unc_y_1,y_e_1 - 1.96*Unc_y_1, alpha=0.5)

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Model_1_Predictions.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 


# Plot quadratic model with Uncertainties
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.scatter(x_s,y_s,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(xg,y_e_2,'r--',linewidth=2)
plt.fill_between(xg, y_e_2 + 1.96*Unc_y_2,y_e_2 - 1.96*Unc_y_2, alpha=0.5)

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Model_2_Predictions.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 



# Plot Cubic model with Uncertainties
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.scatter(x_s,y_s,c='black',
            marker='o',edgecolor='black',s=16)
plt.plot(xg,y_e_3,'r--',linewidth=2)
plt.fill_between(xg, y_e_3 + 1.96*Unc_y_3,y_e_3 - 1.96*Unc_y_3, alpha=0.5)

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)
Name='Model_3_Predictions.png'
plt.tight_layout()
plt.savefig(Name, dpi=200) 

