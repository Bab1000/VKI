# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 10:09:57 2022

@author: mendez
"""


# Adapted from here 
# https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html

import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA


plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

from sklearn import manifold

# READEM: you should run first "Get Cylinder DATA"

#%% Load the data
data = np.load('Snapshot_Matrices.npz')
X_S=data['X_S']; Y_S=data['Y_S']
Xg=data['Xg']; Yg=data['Yg']; t=data['t']; Fs=1/(t[2]-t[1]); n_t=len(t)
D_U=data['D_U']; D_V=data['D_V']
# Assembly both components in one single matrix
D=np.concatenate([D_U,D_V],axis=0) # Reshape and assign


# We look for the decomposition in time:
# in the modal decomposition context, we look for the psi

# The dataset is available at Fs=3000 Hz. We consider a 600 Hz sampling
STEP=5; LIM=10000
X=D[:,0:LIM:STEP].T # We skip some steps 
n_features=np.shape(X)[1]
#% Rescale for convinience
X_scale=(X-X.min())/(X.max()-X.min())
 # Just to track the time.... we use it as a color!
color=t[0:LIM:STEP]


#%% 1.Perform the PCA
pca = PCA(n_components=3)
pca.fit(X_scale)
z = pca.transform(X_scale)

#Check error
X_auto_E=pca.inverse_transform(z)
Err=np.linalg.norm(X_auto_E-X_scale)/np.linalg.norm(X_scale)
print('Error PCA: '+str(Err*100)+' %')


#% 2D plot (PCA)

fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('PCA/POD Manifold in 2D ',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='Cyli_in_Z_2D_PCA.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()


#%% 2. Perform the kPCA
from sklearn.decomposition import KernelPCA
kpca = KernelPCA(kernel="rbf", n_components=3,
                 gamma=0.2,fit_inverse_transform=True,alpha=2,random_state=2)
kpca.fit_transform(X_scale)
z = kpca.transform(X_scale)

#Check the reconstruction error and don't be fooled by the result
X_auto_E=kpca.inverse_transform(z)
Err=np.linalg.norm(X_auto_E-X_scale)/np.linalg.norm(X_scale)
print('Error kPCA: '+str(Err*100)+' %')


#% 3D plot (kPCA)
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(111, projection='3d')
ax.dist = 12
p=ax.scatter(z[0:len(z):2, 0],
           z[0:len(z):2, 1],
           z[0:len(z):2, 2], c=color[1:len(z):2],
           cmap=plt.cm.Set1, edgecolor='k', s=40)
ax.set_xlabel('$\psi_{\mathcal{K}1}$',fontsize=16)
ax.set_ylabel('$\psi_{\mathcal{K}2}$',fontsize=16)
ax.set_zlabel('$\psi_{\mathcal{K}3}$',fontsize=16)
#plt.tight_layout()
Name='kPCA_Manifold_cylinder.png'
plt.colorbar(p)
plt.savefig(Name, dpi=300) 
plt.show()



fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('kPCA Manifold in 2D ',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='Cyli_in_Z_2D_kPCA.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()




#%% 3. Perform the ANN
from keras.layers import Input, Dense
from keras.models import Model
# from keras.layers import  Dropout,BatchNormalization
import tensorflow as tf
tf.random.set_seed(221) ## Remove this to experience the randomness!!!!

from tensorflow import keras
keras.backend.clear_session()
dim=n_features # This is the original dimensionality
encoding_dim = 3  # 2 This is where the compression will take place

# Construct the ANN ()
# this is our input placeholder
input_img = Input(shape=(dim,))
# "encoded" is the encoded representation of the input
encoded = Dense(encoding_dim, activation='linear')(input_img)
#encoded=BatchNormalization(momentum=0.99)(input_img)
#encoded=Dropout(0.2)(encoded)
# "decoded" is the lossy reconstruction of the input
decoded = Dense(dim, activation='linear')(encoded)
# this model maps an input to its reconstruction
autoencoder = Model(input_img, decoded)
# this model maps an input to its encoded representation
encoder = Model(input_img, encoded)

# We compile the autoencoder
autoencoder.compile(optimizer='adam', loss='mse')
# Get a summary
autoencoder.summary()

##%% We do the splitting into test and train sets
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(X_scale, 
                                                    X_scale, 
                                                    test_size=0.1) 
#%% Run the training
History=autoencoder.fit(x_train, x_train,
                epochs=200,
                batch_size=512,
                shuffle=True,
                validation_data=(x_test, x_test))

loss_linlin=History.history['loss']
loss_v=History.history['val_loss']


z=encoder.predict(X_scale)

#Check error
X_auto_E=autoencoder(X_scale)
Err=np.linalg.norm(X_auto_E-X_scale)/np.linalg.norm(X_scale)
print('Error ANN: '+str(Err*100)+' %')



#% 3D plot (ANN)
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(111, projection='3d')
ax.dist = 12
p=ax.scatter(z[0:len(z):2, 0],
           z[0:len(z):2, 1],
           z[0:len(z):2, 2], c=color[1:len(z):2],
           cmap=plt.cm.Set1, edgecolor='k', s=40)
ax.set_xlabel('$\psi_{\mathcal{A}1}$',fontsize=16)
ax.set_ylabel('$\psi_{\mathcal{A}2}$',fontsize=16)
ax.set_zlabel('$\psi_{\mathcal{A}3}$',fontsize=16)
#plt.tight_layout()
Name='ANN_Manifold_cylinder.png'
plt.colorbar(p)

plt.savefig(Name, dpi=300) 
plt.show()




fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('ANN Manifold in 2D ',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='Cyli_in_Z_2D_ANN.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()




#%% 4. Perform the LLE
LLE = manifold.LocallyLinearEmbedding(n_components=3,
                         random_state=56,
                         n_neighbors=15)
z = LLE.fit_transform(X_scale)

# look at the orbits
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(111, projection='3d')
ax.dist = 12
p=ax.scatter(z[0:len(z):2, 0],
           z[0:len(z):2, 1],
           z[0:len(z):2, 2], c=color[1:len(z):2],
           cmap=plt.cm.Set1, edgecolor='k', s=40)
ax.set_xlabel('$\psi_{\mathcal{L}1}$',fontsize=16)
ax.set_ylabel('$\psi_{\mathcal{L}2}$',fontsize=16)
ax.set_zlabel('$\psi_{\mathcal{L}3}$',fontsize=16)
#plt.tight_layout()
Name='LLE_Manifold_cylinder.png'
plt.colorbar(p)
plt.savefig(Name, dpi=300) 
plt.show()



fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('LLE Manifold in 2D ',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='Cyli_in_Z_2D_LLE.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()




#%% 5. Perform the t-SNE in a 2D field
tsne = manifold.TSNE(n_components=2, 
                         perplexity=30,
                         random_state=2,
                         min_grad_norm=1e-06,
                         learning_rate=200,
                         n_iter_without_progress=300,
                         init='pca',
                         early_exaggeration=12,
                         verbose=100)
z = tsne.fit_transform(X_scale)



fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('T-SNE Manifold in 2D ',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='Cyli_in_Z_2D_tSNE.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()


#%% Perform on a 3D plane
tsne = manifold.TSNE(n_components=3, 
                         perplexity=30,
                         random_state=2,
                         min_grad_norm=1e-05,
                         learning_rate=200,
                         n_iter_without_progress=300,
                         init='pca',
                         early_exaggeration=12,
                         verbose=100)
z = tsne.fit_transform(X_scale)


fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(111, projection='3d')
ax.dist = 12
p=ax.scatter(z[0:len(z):2, 0],
           z[0:len(z):2, 1],
           z[0:len(z):2, 2], c=color[1:len(z):2],
           cmap=plt.cm.Set1, edgecolor='k', s=40)
ax.set_xlabel('$\psi_{\mathcal{S}1}$',fontsize=16)
ax.set_ylabel('$\psi_{\mathcal{S}2}$',fontsize=16)
ax.set_zlabel('$\psi_{\mathcal{S}3}$',fontsize=16)
#plt.tight_layout()
Name='tSNE_Manifold_cylinder.png'
plt.colorbar(p)
plt.savefig(Name, dpi=300) 
plt.show()



