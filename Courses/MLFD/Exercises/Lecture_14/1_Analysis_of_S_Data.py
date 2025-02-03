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


#%% Load the the S dataset
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

#%% 1.Perform the PCA
from sklearn.decomposition import PCA
pca = PCA(n_components=2); pca.fit(X_scale)
z = pca.transform(X_scale)
x_app=pca.inverse_transform(z)
Err=np.linalg.norm(x_app-X_scale)/np.linalg.norm(X_scale)
print('Error PCA: '+str(Err*100)+' %')


#% 2D plot (PCA)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
ax.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='S_in_Z_PCA.png'
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()



#%% 2. Perform the kPCA
from sklearn.decomposition import KernelPCA
kpca = KernelPCA(kernel="rbf", n_components=2, 
                 gamma=0.01, fit_inverse_transform=True)
kpca.fit_transform(X_scale)
z = kpca.transform(X_scale)

#Check the reconstruction error after transform
x_app=kpca.inverse_transform(z)
Err=np.linalg.norm(x_app-X_scale)/np.linalg.norm(X_scale)
print('Error kPCA: '+str(Err*100)+' %')


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(x_app[:, 0], x_app[:, 1], x_app[:, 2], c=color, cmap=plt.cm.Spectral)
# ax.view_init(10, -72)
# plt.savefig('S_dataset.png',dpi=300)


#% 2D plot (kPCA)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
ax.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('kPCA in z',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='S_in_Z_kPCA.png'
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()


#%% 3. Perform the ANN
from keras.layers import Input, Dense
from keras.models import Model
#from keras.layers import  Dropout,BatchNormalization
from tensorflow import keras
keras.backend.clear_session()

dim=3 # This is the original dimensionality
encoding_dim = 2  # 2 This is where the compression will take place
# Construct the ANN ()
# this is our input placeholder
input_img = Input(shape=(dim,))
# "encoded" is the encoded representation of the input
encoded = Dense(encoding_dim, activation='elu')(input_img)
# "decoded" is the lossy reconstruction of the input
decoded = Dense(dim, activation='elu')(encoded)
# this model maps an input to its reconstruction
autoencoder = Model(input_img, decoded)
# this model maps an input to its encoded representation
encoder = Model(input_img, encoded)

# We compile the autoencoder
autoencoder.compile(optimizer='adam', loss='mse')
# Get a summary
autoencoder.summary()
##% We do the splitting into test and train sets
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(X_scale,X_scale, 
                                                    test_size=0.1) 
#% Run the training
History=autoencoder.fit(x_train, x_train,
                epochs=1000,
                batch_size=128,
                shuffle=True,
                validation_data=(x_test, x_test))
# get convergence history
loss_linlin=History.history['loss']
loss_v=History.history['val_loss']

# Prediction of the encoding
z=encoder.predict(X_scale)
x_tilde=autoencoder.predict(X_scale)
#Check error
Err=np.linalg.norm(x_tilde-X_scale)/np.linalg.norm(X_scale)
print('Error ANN: '+str(Err*100)+' %')





#% 2D plot (ANN)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
ax.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('ANN Autoencoder in z',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='S_in_Z_ANN.png'
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()




#%% 4. Perform the LLE
from sklearn import manifold
LLE = manifold.LocallyLinearEmbedding(n_components=2,
                         random_state=56,
                         n_neighbors=30)
z = LLE.fit_transform(X_scale)


#% 2D plot (LLE)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
ax.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('LLE Manifold',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='S_in_Z_LLE.png'
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()




#%% 5. Perform the t-SNE
tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=56,
                         perplexity=50)
z = tsne.fit_transform(X)


#% 2D plot (t_sne)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
ax.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('t-SNE Manifold in z',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='S_in_Z_tSNE.png'
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()








