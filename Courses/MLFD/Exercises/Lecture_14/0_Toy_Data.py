# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:08:41 2024

@author: mendez
"""


import numpy as np
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

#%% Load the S Dataset
from sklearn import datasets

n_points = 2000
X, color = datasets.make_s_curve(n_points, random_state=0)
X_scale=(X-X.min())/(X.max()-X.min())
# Visualize the data
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X_scale[:, 0], X_scale[:, 1], X_scale[:, 2], c=color, cmap=plt.cm.Spectral)
ax.view_init(10, -72)
plt.savefig('S_dataset.png',dpi=300)

#%% Load the MNIST dataset
digits = datasets.load_digits(n_class=6)
X, y = digits.data, digits.target
n_samples, n_features = X.shape
# We will attribute the variable y to 'color'
color=y

#% Plot the first 100 digits
fig, axs = plt.subplots(nrows=10, ncols=10, figsize=(6, 6))
for idx, ax in enumerate(axs.ravel()):
    ax.imshow(X[idx].reshape((8, 8)), cmap=plt.cm.binary)
    ax.axis("off")
_ = fig.suptitle("A selection from the 64-dimensional digits dataset", fontsize=16)
plt.savefig('MNIST.png', dpi=300) 
plt.show()


#%% Check one of them... 8x8 is not a good resolution... is it :D?
fig = plt.figure(figsize=(4,4))
Index=74
plt.imshow(X[Index,:].reshape(8,8),cmap=plt.cm.binary)
plt.savefig('MIST_'+str(Index)+'.png', dpi=300) 
plt.show()

fig = plt.figure(figsize=(4,4))
Index=33
plt.imshow(X[Index,:].reshape(8,8),cmap=plt.cm.binary)
plt.savefig('MIST_'+str(Index)+'.png', dpi=300) 
plt.show()


fig = plt.figure(figsize=(4,4))
Index=521
plt.imshow(X[Index,:].reshape(8,8),cmap=plt.cm.binary)
plt.savefig('MIST_'+str(Index)+'.png', dpi=300) 
plt.show()

#% Rescale for convinience
X_scale=(X-X.min())/(X.max()-X.min())