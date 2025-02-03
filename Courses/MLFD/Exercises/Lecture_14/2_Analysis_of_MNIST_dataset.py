# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:08:41 2024

@author: mendez
"""


# Adapted from here 
# https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html

import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA, KernelPCA

plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

from sklearn import manifold

#%% Load the data
from sklearn import datasets
digits = datasets.load_digits(n_class=6)
X, y = digits.data, digits.target
n_samples, n_features = X.shape
# We will attribute the variable y to 'color'
color=y

#% Rescale for convinience
X_scale=(X-X.min())/(X.max()-X.min())
 


#%% 1.Perform the PCA--------------------------
n_neighbors = 30
pca = PCA(n_components=2)
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
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='MNIST_in_Z_PCA.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()



#%% 2. Perform the kPCA with various gammas 
Gammas=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7])

for k in range(len(Gammas)):
    kpca = KernelPCA(kernel="rbf", n_components=2,
                 gamma=Gammas[k],fit_inverse_transform=True,alpha=1)
    kpca.fit_transform(X_scale)
    z = kpca.transform(X_scale)
    #Check the reconstruction error after transform

    #Check error
    X_auto_E=kpca.inverse_transform(z)
    Err=np.linalg.norm(X_auto_E-X_scale)/np.linalg.norm(X_scale)
    print('Error kPCA: '+str(Err*100)+' %')

    #% 2D plot (kPCA)
    fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
    plt.scatter(z [:, 0], z [:, 1], c=color,
               cmap=plt.cm.Set1, edgecolor='k', s=40)
    plt.title('kPCA in z',fontsize=14)
    ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
    ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
    Name='MNIST_in_Z_kPCA'+str(Gammas[k])+'.png'
    plt.colorbar()
    
    plt.tight_layout()
    plt.savefig(Name, dpi=300) 
    plt.show()



#%% 3. Perform the ANN
from keras.layers import Input, Dense
from keras.models import Model
# from keras.layers import  Dropout,BatchNormalization

from tensorflow import keras
keras.backend.clear_session()
dim=n_features # This is the original dimensionality
encoding_dim = 2  # 2 This is where the compression will take place

# Construct the ANN ()
# this is our input placeholder
input_img = Input(shape=(dim,))
# "encoded" is the encoded representation of the input
encoded = Dense(encoding_dim, activation='elu')(input_img)
#encoded=BatchNormalization(momentum=0.99)(encoded)
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

##%% We do the splitting into test and train sets
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(X_scale, 
                                                    X_scale, 
                                                    test_size=0.1) 
#%% Run the training
History=autoencoder.fit(x_train, x_train,
                epochs=1000,
                batch_size=128,
                shuffle=True,
                validation_data=(x_test, x_test))

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
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)

plt.title('ANN Autoencoder in z',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='MNIST_in_Z_ANN.png'
plt.colorbar()

plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()



#%% 4. Perform the LLE
LLE = manifold.LocallyLinearEmbedding(n_components=2,
                         random_state=5,
                         n_neighbors=20)
z = LLE.fit_transform(X)


#% 2D plot (LLE)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('LLE Manifold',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='MNIST_in_Z_LLE.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()


#%% 5. Perform the t-SNE
tsne = manifold.TSNE(n_components=2, 
                         perplexity=30,
                         learning_rate=200,
                         n_iter_without_progress=300,
                         init='random',
                         early_exaggeration=12,
                         verbose=100)
z = tsne.fit_transform(X_scale)



#% 2D plot (t_sne)
fig, ax = plt.subplots(figsize=(6, 4)) # This creates the figure
plt.scatter(z [:, 0], z [:, 1], c=color,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
plt.title('t-SNE Manifold in z',fontsize=14)
ax.set_xlabel('$\mathbf{z}[1]$',fontsize=16)
ax.set_ylabel('$\mathbf{z}[2]$',fontsize=16)
Name='MNIST_in_Z_tSNE.png'
plt.colorbar()
plt.tight_layout()
plt.savefig(Name, dpi=300) 
plt.show()








