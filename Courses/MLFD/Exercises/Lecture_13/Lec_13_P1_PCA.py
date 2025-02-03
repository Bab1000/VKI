"""
Created on Sun May 26 09:51:32 2024

@author: mendez
"""

"""
# Tutorial on PCA and Eigenfaces

This tutorial is adapted from https://scipy-lectures.org/packages/scikit-learn/auto_examples/plot_eigenfaces.html

We perform the PCA "manually" to ensure the full understanding of the material presented in class. 
We use both the "classic PCA" and the "dual PCA".

We compute the principal component of a dataset collecting pictures of Olivetti's employees.
These will be eigenfaces. 
The dataset is a compressed version of the original one, but still large 
enough to illustrate the benefits of dimensionality reduction
"""

# Import the necessary libraries and download the dataset
from sklearn import datasets
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
import numpy as np
from scipy.sparse.linalg import svds, eigsh

# This is a customization for the plots
plt.rc('text', usetex=False)      # This is Miguel's customization
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=18)
plt.rc('ytick',labelsize=18)

faces = datasets.fetch_olivetti_faces()
faces.data.shape

# plot some of the faces
fig = plt.figure(figsize=(8, 6))
# plot several images
for i in range(15):
    ax = fig.add_subplot(3, 5, i + 1, xticks=[], yticks=[])
    ax.imshow(faces.images[i], cmap=plt.cm.bone)

# Explore the dataset. This is the image of employ number 0:
image=faces.images[0]
fig=plt.figure(figsize=(4,4))
plt.imshow(image,cmap=plt.cm.bone)
plt.gca().set_aspect('equal')

#%% Split into training and testing data
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(faces.data,
        faces.target, random_state=0)
print(X_train.shape, X_test.shape)

# To use the formulation in the slides, we will transpose these !
X_train=X_train.T
X_test=X_test.T

# To visualize one of these we need to reshape the vector.
# For example, this is the 6th picture in the training set:
image=X_train[:,5].reshape(64,64)
fig=plt.figure(figsize=(4,4))
plt.imshow(image,cmap=plt.cm.bone)
plt.gca().set_aspect('equal')

"""

 Compute the Eigenfaces on the training data

The standard PCA would generally require mean subtraction before decomposition.
We ignore this for the moment (try doing it as an exercise... especially if you
                               compare with the scikit learn implementation!)

We need to solve the following eigenvalue problem:

$$(𝐗 𝐗^T)\tilde{𝐔}=\tilde{𝐔} \tilde{𝚲}\rightarrow 𝐂 \tilde{𝐔}=\tilde{𝐔} \tilde{𝚲} $$


"""

# Let compute the first 150 principal components
n_C=150
# We first compute C:
C=X_train@X_train.T
# Eigenvalue solver for Hermitian Matrices
Lambda, U = eigsh(C, k=n_C)
# It turns out that this does not rank them in decreasing order.
# Hence we do it manually:
idx = np.flip(np.argsort(Lambda)); Lambda = Lambda[idx] ; U = U[:, idx]

# Plot the eigenvalues, measuring the relative contribution of each component
fig, ax = plt.subplots(figsize=(7,6))
plt.plot(Lambda/4096,'ko')
plt.xlabel('$r $',fontsize=22)
plt.ylabel('$\lambda_r/n_s $',fontsize=22)
plt.yscale("log")
plt.tight_layout()
NameOUT='Exercise_1_Lambdas.png'
plt.savefig(NameOUT, dpi=100)
plt.show()


fig = plt.figure(figsize=(8, 6))
for i in range(15):
    ax = fig.add_subplot(3, 5, i + 1, xticks=[], yticks=[])
    ax.imshow(U[:,i].reshape(64,64), cmap=plt.cm.bone)
plt.tight_layout()
NameOUT='Exercise_2_U_s.png'
plt.savefig(NameOUT, dpi=100)
plt.show()

"""# Autoencoding for with different ranks
We build an approximation with 10, 50 , 150 eigenfaces!

"""

R_CUT=50 # The rank at which we cut. Run this with 10, 50, 150
# Pick a random element of the test data
Sample=X_train[:,44]
# Define the reduced set
U_tilde=U[:,0:R_CUT-1]

# This would be the encoding
z=U_tilde.T.dot(Sample)

# Show the encoding
fig, ax = plt.subplots(figsize=(7,4))
plt.stem(z,'ko')
plt.tight_layout()
NameOUT='Exercise_1_Encoding.png'
plt.savefig(NameOUT, dpi=100)
plt.show()

# This would be the decoder
x_tilde=U_tilde.dot(z)

# This would be the autoencoder in one shot:
#A=U_tilde@U_tilde.T
#x_tilde_A=A.dot(Sample)

# Show the Original
image=Sample.reshape(64,64)


fig=plt.figure(figsize=(4,4))
plt.imshow(image,cmap=plt.cm.bone)
plt.gca().set_aspect('equal')
plt.tight_layout()
NameOUT='Exercise_1_Sample.png'
plt.savefig(NameOUT, dpi=100)
plt.show()

# Show the Approximation
image=x_tilde.reshape(64,64)
fig=plt.figure(figsize=(4,4))
plt.imshow(image,cmap=plt.cm.bone)
plt.gca().set_aspect('equal')
plt.tight_layout()
NameOUT='Exercise_1_Approximation_'+str(R_CUT)+'.png'
plt.savefig(NameOUT, dpi=100)
plt.show()


#

""" Dual PCA

We compute the encoding Z from the dual PCA.

It is left as an exercise to show what would be the autoencoding in the dual formalism!


"""

# We now compute the encoding
K=X_train.T@X_train
# Eigenvalue solver for Hermitian Matrices
Lambda_V, V = eigsh(K, k=n_C)
# It turns out that this does not rank them in decreasing order.
# Hence we do it manually:
idx = np.flip(np.argsort(Lambda_V)); Lambda_V = Lambda_V[idx] ; V = V[:, idx]


# Compare t he Z's for the dual for one of the samples
idx=50

Sample=X_train[:,idx]

z_U=U.T.dot(Sample)
# CAREFUL: There can be a sign indeterminacy!
z_V=np.sqrt(np.diag(Lambda)).dot(V[idx,:].T)

# CAREFUL: the sign of the eigenvectors is... random! so we shoul
# make sure they have the same sign if we want to use the decoder U/
for j in range(len(z_V)):
  if z_U[j] * z_V[j] > 0:
        z_V[j]=z_V[j]
  else:
        z_V[j]=z_V[j]*-1


# Decoder
x_approx_U=U.dot(z_U).reshape(64,64); x_approx_V=U.dot(z_V).reshape(64,64)

# Create a figure to hold the plots
fig, ax = plt.subplots(1, 2)

# Display the first image
ax[0].imshow(x_approx_U,cmap=plt.cm.bone)
ax[0].axis('off')  # Turn off axis numbering and ticks
ax[0].set_title('Decoding from z_U')  # Set a title, if desired

# Display the second image
ax[1].imshow(x_approx_V,cmap=plt.cm.bone)
ax[1].axis('off')
ax[1].set_title('Decoding from z_V')

# Display the figure
plt.show()

np.shape(V[idx,:].T)

Check=55
print(z_U[Check],z_V[Check])

