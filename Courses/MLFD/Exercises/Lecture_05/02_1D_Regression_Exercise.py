# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 11:52:32 2022

@author: mendez, ratz
"""

# Basics needed to create and plot the data
import numpy as np
import matplotlib.pyplot as plt
import os

# For the animation
import imageio

# Needed for the feature scaling
from sklearn.preprocessing import MinMaxScaler
# Import the in-built sklearn functions for regression
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

# Configuration for plots  
fontsize = 12
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'serif'
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['legend.fontsize'] = fontsize-3
plt.rcParams['font.size'] = fontsize

# Define the output folder
Fol_Out = '02_1D_Regression_Exercise'
if not os.path.exists(Fol_Out):
    os.makedirs(Fol_Out)

#%% Define all of the bases and the functions to assemble the basis matrices

# The Gaussian
def Gauss_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r
    phi_r = np.exp(-c_r**2*d**2)
    return phi_r

# The Matern C4 kernel
def Matern_C4_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = np.abs(x-x_r)
    phi_r = np.exp(-c_r*d) * (c_r**2*d**2 + 3*c_r*d + 3)
    return phi_r

# The (custom) C4
def Compact_C4_RBF(x, x_r=0, c_r=0.1):
    # Get distance
    d = x-x_r 
    phi_r = (1+d/c_r)**5 * (1-d/c_r)**5
    phi_r[np.abs(d)>c_r] = 0
    return phi_r


#%% Generate and plot the data

# X-values with a gap in between
x1 = np.linspace(0, 4.3, 200, endpoint=True)
x2 = np.linspace(4.8, 10, 200, endpoint=True)
x = np.concatenate((x1, x2))
# Create Y-values with noise
y_clean = 3*x + (x/100)**3 + 4*np.sin(3/2*np.pi*x)
np.random.seed(0)
y = y_clean + 1*np.random.randn(len(x))
# Introduce some outliers in x=2 and x=8
G1 = 10 * np.exp(-(x-2)**2/0.005) * np.random.randn(len(x))
G2 = 15 * np.exp(-(x-8)**2/0.005) * np.random.randn(len(x))
y_final = y + G1 + G2

# Plot the data
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
ax.scatter(x, y_final, c='black', marker='o', edgecolor='black', s=16)
ax.set_xlabel('x')
ax.set_ylabel('y')
Name = 'Exercise_1_data.png'
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + Name, dpi=300) 

plt.show()

#%% Perform the feature scaling

# Scale the x-values
scaler_X = MinMaxScaler() 
scaler_X.fit_transform(x.reshape(-1, 1))
x_prime = scaler_X.transform(x.reshape(-1, 1))

# Scale also the y-values
scaler_Y = MinMaxScaler() 
scaler_Y.fit_transform(y_final.reshape(-1, 1))
y_prime = scaler_Y.transform(y_final.reshape(-1, 1))


# Plot the data
fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
ax.scatter(x_prime, y_prime, c='black', marker='o', edgecolor='black', s=16)
ax.set_xlabel('x')
ax.set_ylabel('y')
Name = 'Exercise_1_data_scaled.png'
fig.tight_layout()
fig.savefig(Fol_Out + os.sep + Name, dpi=300) 

#%% These are the functions to assemble the entire basis matrix Phi

# The Gaussian
def PHI_Gauss_X(x_in, x_b, c_r=0.05):
    n_x = np.size(x_in)
    n_b = np.size(x_b)
    # Initialize Basis Matrix on x
    Phi_X = np.zeros((n_x, n_b+1))
    # Add a linear term
    Phi_X[:, 0] = x_in 
    # Fill the rest of the matrix with Gaussians
    for j in range(0, n_b): 
        Phi_X[:, j+1] = Gauss_RBF(x_in, x_r=x_b[j], c_r=c_r)  
    return Phi_X

# The Matern C4
def PHI_Matern_C4_X(x_in, x_b, c_r=0.05):
    n_x = np.size(x_in)
    n_b = np.size(x_b)
    # Initialize Basis Matrix on x
    Phi_X = np.zeros((n_x, n_b+1))
    # Add a linear term
    Phi_X[:, 0] = x_in 
    # Fill the rest of the matrix with sigmoids
    for j in range(0, n_b):
        Phi_X[:, j+1] = Matern_C4_RBF(x_in, x_r=x_b[j], c_r=c_r)
    return Phi_X

# The (custom) C4
def PHI_C4_X(x_in, x_b, c_r=0.1):
    n_x = np.size(x_in)
    n_b = len(x_b)
    # Initialize Basis Matrix on x
    Phi_X = np.zeros((n_x, n_b+1))
    # Add a linear term
    Phi_X[:, 0] = x_in 
    # Fill the rest of the matrix with Gaussians
    for j in range(0, n_b):
        Phi_X[:, j+1] = Compact_C4_RBF(x_in, x_r=x_b[j], c_r=c_r)
    return Phi_X

#%% Question 1

# Define the number of bases
n_b = 100
# Define grid of collocation points
x_b = np.linspace(0, 1, n_b)
    

# Check the Gaussian basis:
Phi_Gauss_X = PHI_Gauss_X(x_prime[:,0], x_b, c_r=1/0.05)
plt.figure(dpi=200)
plt.plot(x_prime, Phi_Gauss_X[:, ::10])

# Compute the condition number 
H_Gauss = Phi_Gauss_X.T@Phi_Gauss_X
k_H_Gauss = np.linalg.cond(H_Gauss)
print('rcond for Gauss RBF: {:.3f}'.format(k_H_Gauss))  


# Check the Sigmoid basis:
PHI_MC4_X = PHI_Matern_C4_X(x_prime[:, 0], x_b, c_r=60)
plt.figure(dpi=200)
plt.plot(x_prime, PHI_MC4_X[:, ::10])

# Compute the condition number 
H_MC4 = PHI_MC4_X.T@PHI_MC4_X
k_H_MC4 = np.linalg.cond(H_MC4)
print('rcond for Matern C4: {:.3f}'.format(k_H_MC4))      


# Check the C4 basis:
Phi_C4_X = PHI_C4_X(x_prime[:,0], x_b, c_r=0.1)
plt.figure(dpi=200)
plt.plot(x_prime, Phi_C4_X[:, ::10])

# Compute the condition number 
H_C4 = Phi_C4_X.T@Phi_C4_X
k_H_C4 = np.linalg.cond(H_C4)
print('rcond for C4 RBF: {:.3f}'.format(k_H_C4))   

plt.show() 

#%% Question 2

# Ordinary Least Square (OLS) regression using numpy's inversion
H_C4 = Phi_C4_X.T@Phi_C4_X
w_C4 = np.linalg.inv(H_C4).dot(Phi_C4_X.T).dot(y_prime)   

# Make prediction on the nex x (regularly spaced)
x_test = np.linspace(0, 1, 300)
Phi_X_test = PHI_C4_X(x_test, x_b, c_r=0.1)
y_prime_pred = Phi_X_test.dot(w_C4) 

# Alternative: Use sklearn's inbuilt function
reg = LinearRegression(fit_intercept=False).fit(Phi_C4_X, y_prime)
w_C4_scikit = reg.coef_.T  # Weights can be extracted for comparison
y_prime_pred_scikit = reg.predict(Phi_X_test)


# Plot the numpy and scikit regressions on a global scale
fig, ax = plt.subplots(figsize=(5, 3), dpi=200) 
ax.scatter(x_prime, y_prime, c='white', marker='o', edgecolor='black',
            s=10, label='Data')
ax.plot(x_test, y_prime_pred, 'r--', label='numpy')
ax.plot(x_test, y_prime_pred_scikit, 'b--', label='scikit')
ax.set_xlabel('x')
ax.set_ylabel('y') 
ax.legend()
fig.tight_layout()
Name = 'NPvsScikit_OLS.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

# Plot the numpy and scikit regressions on a zoomed scale
fig, ax = plt.subplots(figsize=(5, 3), dpi=200) 
ax.scatter(x_prime, y_prime, c='white', marker='o', edgecolor='black',
           s=10, label='Data')
ax.plot(x_test, y_prime_pred, 'r--', label='numpy')
ax.plot(x_test, y_prime_pred_scikit, 'b--', label='scikit')
ax.set_xlabel('x')
ax.set_ylabel('y')  
ax.legend()
ax.set_xlim([0.5, 0.7])
ax.set_ylim([0, 1])
fig.tight_layout()
Name = 'NPvsScikit_OLS_zoom.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

plt.show()

#%% Question 3 

# Ordinary Least Square
OLS = LinearRegression(fit_intercept=False) 
OLS.fit(Phi_C4_X, y_prime)
w_reg_s = OLS.coef_
y_prime_pred_O = OLS.predict(Phi_X_test)

# Lasso Regression
Las = Lasso(fit_intercept=False, alpha=0.01)
Las.fit(Phi_C4_X, y_prime)
w_reg_L_s = Las.coef_
y_prime_pred_L = Las.predict(Phi_X_test)

# Ridge Regression
Rid = Ridge(fit_intercept=False, alpha=0.01) 
Rid.fit(Phi_C4_X, y_prime)
w_reg_R_s = Rid.coef_
y_prime_pred_R = Rid.predict(Phi_X_test)

# Plot the three on a global scale
fig, ax = plt.subplots(figsize=(5, 3), dpi=200) 
ax.scatter(x_prime, y_prime, c='white', marker='o', edgecolor='black',
            s=10, label='Data')
ax.plot(x_test, y_prime_pred_O, 'r--', label='OLS')
ax.plot(x_test, y_prime_pred_R, 'b--', label='Ridge')
ax.plot(x_test, y_prime_pred_L, 'k--', label='Lasso')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend()
fig.tight_layout()
Name = 'Lasso_Ridge_OLS_full.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

# Plot the three on a zommed scale
fig, ax = plt.subplots(figsize=(5, 3), dpi=200) 
ax.scatter(x_prime,y_prime,c='white',
            marker='o',edgecolor='black',
            s=10,label='Data')
ax.plot(x_test,y_prime_pred_O,'r--',label='OLS')
ax.plot(x_test,y_prime_pred_R,'b--',label='Ridge')
ax.plot(x_test,y_prime_pred_L,'k--',label='Lasso')
ax.set_xlabel('x')
ax.set_ylabel('y')  
ax.legend()
ax.set_xlim([0.5, 0.7])
ax.set_ylim([0, 1])
fig.tight_layout()
Name = 'Lasso_Ridge_OLS_Zoom.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

plt.show()

#%% Plot the weights

# OLS
fig, ax = plt.subplots(figsize=(5, 2), dpi=200) 
ax.set_title('Weights from OLS',fontsize=16)
ax.stem(w_reg_s[0, :])
ax.set_xlabel('$r$')
ax.set_ylabel('$w_r$')  
fig.tight_layout()
Name = 'W_OLS.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

# Ridge
fig, ax = plt.subplots(figsize=(5, 2), dpi=200) 
ax.set_title('Weights from Ridge R',fontsize=16)
ax.stem(w_reg_R_s[0, :])
ax.set_xlabel('$r$')
ax.set_ylabel('$w_r$')
fig.tight_layout()
Name = 'W_R.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

# OLS
fig, ax = plt.subplots(figsize=(5, 2), dpi=200) 
ax.set_title('Weights from Lasso R')
ax.stem(w_reg_L_s)
ax.set_xlabel('$r$')
ax.set_ylabel('$w_r$')  
fig.tight_layout()
Name = 'W_L.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200) 

plt.show()

#%% Exercise 4

# Here we perform the K-fold validation using sklearn
K = 5 
kf = KFold(n_splits=K, random_state=3, shuffle=True)

# counter to save the folds
count = 1
# Have a look at the data spliting!
for train_index, test_index in kf.split(Phi_C4_X):
    Phi_X_train, Phi_X_test = Phi_C4_X[train_index], Phi_C4_X[test_index]
    fig, ax = plt.subplots(figsize=(5, 3), dpi=200)  
    y_train, y_test = y_prime[train_index], y_prime[test_index]
    x_train, x_test = x_prime[train_index], x_prime[test_index]
    ax.plot(x_train, y_train, 'bo', label='Train')
    ax.plot(x_test, y_test, 'ro', label='Test')
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('y')  
    fig.tight_layout()
    Name = 'Fold_' + str(count) + '.png'
    count+=1
    fig.savefig(Fol_Out + os.sep + Name, dpi=200) 


Animation_Name = Fol_Out + os.sep + 'K_Folds.gif'

# Make the animation
images = []
for k in range(1, 5):
    MEX = 'Preparing Im ' + str(k) 
    print(MEX)  
    FIG_NAME = Fol_Out + os.sep + 'Fold_' + str(k) + '.png'  
    images.append(imageio.imread(FIG_NAME))
    
# Assembly the video
imageio.mimsave(Animation_Name, images, duration=0.3)

plt.show()

#%% Question 4

# Here we perform the K-fold validation
alphas = np.linspace(0, 0.001, 501)
J_out = np.zeros(len(alphas))
J_in = np.zeros(len(alphas))
# Create the k-fold split object
kf = KFold(n_splits=K, random_state=3, shuffle=True)

# Loop over given alphas
for j in range(len(alphas)):
    print('alpha ' + str(j) + ' of ' + str(len(alphas)))
    # Select one alpha
    alpha = alphas[j]  
    # Initialize the out of sample error vector
    count = 0
    J_out_fold = np.zeros(K)
    J_in_fold = np.zeros(K)
    # Loop over the folds
    for train_index, test_index in kf.split(Phi_C4_X):
        # Get the training and test sets
        Phi_X_train, Phi_X_test = Phi_C4_X[train_index], Phi_C4_X[test_index] 
        y_train, y_test = y_prime[train_index], y_prime[test_index] 
        # Fit the model on the trainig set
        Rid = Ridge(fit_intercept=False, alpha=alpha) 
        Rid.fit(Phi_X_train, y_train)
        # Test the model on the testing set
        y_prime_train = Rid.predict(Phi_X_train) 
        y_prime_test = Rid.predict(Phi_X_test) 
        # Collect all the out of sample errors
        J_in_fold[count] = np.mean((y_prime_train-y_train)**2) 
        J_out_fold[count] = np.mean((y_prime_test-y_test)**2) 
    # Take the mean out of sample error over the folds
    J_in[j] = np.mean(J_in_fold)
    J_out[j] = np.mean(J_out_fold)
          
fig, ax = plt.subplots(figsize=(5, 3), dpi=200) 
ax.plot(alphas, J_in, label='$J_i$')
ax.plot(alphas, J_out, label='$J_o$')
ax.legend()
ax.set_xlabel('$\\alpha$')
ax.set_ylabel('$J_i, J_o$')  
fig.tight_layout()
Name = 'J_i_J_o.png'
fig.savefig(Fol_Out + os.sep + Name, dpi=200)

plt.show()
