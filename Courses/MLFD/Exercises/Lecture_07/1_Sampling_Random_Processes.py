# -*- coding: utf-8 -*-
"""
Created on Sat May 25 20:24:56 2024

@author: mendez
"""

import numpy as np
import matplotlib.pyplot as plt


# Show the result of the fit 
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)



#%% Kernel Definition
from sklearn.metrics.pairwise import euclidean_distances
def kernel_N(X1,X2,l_c=1,sigma_f=1):
    sqdist = euclidean_distances(X1, X2)**2
    return sigma_f**2 * np.exp(-0.5 / l_c**2 * sqdist)


#%% Krasser's function for plotting Gaussian Processes
def plot_gp(mu, cov, X, X_train=None, Y_train=None, samples=[]):
    #%P=
    X = X.ravel()
    mu = mu.ravel()
    uncertainty = 1.96 * np.sqrt(np.diag(cov))
    
    plt.fill_between(X, mu + uncertainty, mu - uncertainty, alpha=0.1)
    plt.plot(X, mu, label='Mean')
    for i, sample in enumerate(samples):
        plt.plot(X, sample, lw=1, ls='--', label=f'Sample {i+1}')
    if X_train is not None:
        plt.plot(X_train, Y_train, 'rx')
    plt.legend()


#%% Sample grid (given as a vector)
X = np.linspace(0,5,200).reshape(-1, 1)
# Parameters of the process
l_c=0.5; sigma_f=1
# Prepare first and second order statistics (custom functions)
mu = np.zeros(X.shape)
K = kernel_N(X, X,l_c=l_c)
# Scikit implementation
from sklearn.metrics.pairwise import rbf_kernel
K=sigma_f**2*rbf_kernel(X,X,gamma=0.5 / l_c**2)

# Sample the process 3 times and plot the results
markers=['k','r','b']

fig, ax = plt.subplots(figsize=(7, 3)) # This creates the figure
for i in range(3):
  f_x_i = np.random.multivariate_normal(mu.ravel(), K, 1)
  plt.plot(X,f_x_i.T,markers[i],label='sample '+str(i+1))

plt.title('Samples with $l_c=$'+str(l_c),size=18)
plt.legend()
ax.set_xlabel('x',fontsize=18)
ax.set_ylabel('y',fontsize=18)
Name_FIG='3_Samples_l_c'+str(l_c)+'.png'
plt.savefig(Name_FIG)
plt.show()




#%% More sophisticated plot with uncertainties

# Draw several samples
n_p=6
samples = np.random.multivariate_normal(mu.ravel(), K, n_p)

fig, ax = plt.subplots(figsize=(5, 3)) # This creates the figure
plot_gp(mu, K, X, samples=samples)




