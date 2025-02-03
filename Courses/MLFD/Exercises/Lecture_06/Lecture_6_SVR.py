# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:51:05 2025

@author: mendez
"""

# Import essential things
import numpy as np
import matplotlib.pyplot as plt

# plot customization (a matter of taste)

# Configuration for plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)



#%% Section 0: Create the usual dataset
x1 = np.linspace(0, 4.3, 200, endpoint=True)
x2 = np.linspace(4.8, 10, 200, endpoint=True)
x_data=np.concatenate((x1,x2))
# Create the deterministic part
y_clean= 3*x_data+(x_data/100)**3+4*np.sin(3/2*np.pi*x_data)
# Add (a seeded) stochastic part
np.random.seed(0)
y=y_clean+1*np.random.randn(len(x_data))
# Introduce some outliers in x=2 and x=8
G1=10*np.exp(-(x_data-2)**2/0.005)*np.random.randn(len(x_data))
G2=15*np.exp(-(x_data-8)**2/0.005)*np.random.randn(len(x_data))
y_data=y+G1+G2


#%% Support vector regression in python
from sklearn.svm import SVR
# New predictions requested here:
x_ss=np.linspace(x_data.min(),x_data.max(),500)
# Create SVR object
svr=SVR(kernel='rbf',gamma=10,C=50,epsilon=2)
# Fit the regressor
svr.fit(x_data.reshape(-1,1),np.ravel(y_data))
# Make predictions:
y_p_SVM=svr.predict(x_ss.reshape(-1,1))    
# Look for the epsilon sensitive tube
y_p_eps=y_p_SVM+svr.epsilon
y_m_eps=y_p_SVM-svr.epsilon

# Plot the result
fig, ax = plt.subplots(figsize=(5, 3)) 
plt.scatter(x_data,y_data,c='white',marker='o',edgecolor='black',
            s=10,label='Data')
plt.plot(x_ss,y_p_SVM,'r--')
plt.plot(x_ss,y_m_eps,'b--')
plt.plot(x_ss,y_p_eps,'b--')

ax.set_xlabel('x',fontsize=16)
ax.set_ylabel('y',fontsize=16)  

Name='SVR.png'

plt.tight_layout()
plt.savefig(Name, dpi=200) 







