# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 16:32:04 2023

@author: ahizi
"""

# import torch
import matplotlib.pyplot as plt
import numpy as np
import time
from scipy.integrate import solve_ivp


#%%

# lorentz = np.float32(np.expand_dims(np.load('states_lorentz.npy').T,axis=1))

def Lorentz(t, y, sigma=10.0, rho=28.0, beta=8.0/3.0):
        "Lorenz attractor - default coefficients lead to chaotic behaviour"
        x, y, z = y[0], y[1], y[2]
        dxdt = sigma * (y - x)
        dydt = x * (rho - z) - y
        dzdt = x * y - beta * z
        
        return np.hstack((dxdt, dydt, dzdt))

def Spiral(t, y, a=-1/10, b=2, c=-1/10):
        "Spiral - default coefficients lead to spiral converging to origin"
        
        A = np.array([[a, b, 0], [-b, a, 0],[0,0,c]])
        dydt = np.matmul(A,y**3)
        
        return dydt
    
def generateData(y0 = [-5. , -5., 40.] ,dataSize = 1000,t0=0, tf = 15.0, name='lorentz'):
    if name=='lorentz':
        func=Lorentz
    if name=='spiral':
        func=Spiral
    t = np.linspace(t0, tf, dataSize)
    states = solve_ivp(func, (t[0],t[-1]), y0, t_eval=t).y
    np.save(f't_{name}.npy', t)
    np.save(f'states_{name}.npy', states)

def visualize(t, states, prediction=None, iter=0,save=False,gif=False,fz=10):
    fig = plt.figure()
    subfigs = fig.subfigures(1, 2)
    ax1 = subfigs[1].add_subplot(projection='3d')
    ax0 = subfigs[0].add_subplot()
    
    ax1.set_title('3D view',fontsize=fz)
    ax0.set_title('Trajectory',fontsize=fz)
  
    ax1.axes.set_xlim3d(min(states[0,:])*0.8-2, max(states[0,:])*1.2+2) 
    ax1.axes.set_ylim3d(min(states[1,:])*0.8-2, max(states[1,:])*1.2+2)
    ax1.axes.set_zlim3d(min(states[2,:])*0.8-2, max(states[2,:])*1.2+2)
    
    ax1.plot(states[0,:], states[1,:], states[2,:], 'green', label='true trajectory')
    ax0.plot(t,states[0,:], label='true x')
    ax0.plot(t,states[1,:], label='true y')
    ax0.set_ylim(np.amin(states[:2,:])*0.8-2,np.amax(states[:2,:])*1.2+2)
    ax0.set_xlabel('Time',fontsize=fz)
    ax0.set_ylabel('Position',fontsize=fz)
    
    if states[2,0]!=np.mean(states[2,:]):
        ax0.plot(t,states[2,:], label='true z')
        ax0.set_ylim(np.amin(states)*0.8-2,np.amax(states)*1.2+2)
    if prediction is not None:
        
        ax1.scatter(states[0,:], states[1,:], states[2,:],color='black', label='sampled data',s=2)

        
        pred = prediction

        ax1.scatter(pred[0,:], pred[1,:], pred[2,:], color='red',label='Prediction',s=2)
        ax0.plot(t,pred[0,:],':', color= 'C0', label='Prediction')
        ax0.plot(t,pred[1,:],':', color= 'C1')
        ax0.plot(t,pred[2,:],':', color= 'C2')
        ax0.set_title(f'Trajectory, iteration {iter}')
    
       
    ax0.legend(bbox_to_anchor=(1.1, 1.05))
    if save!= False:
        plt.savefig(save,dpi=250)
    plt.show()
    

def make_gif(frame_folder='./plots/', name='RNN'):
    import glob
    import os
    from PIL import Image

    files = list(glob.glob(f"{frame_folder}/{name}*"))
    files.sort(key=os.path.getctime)
    print('Saving gif...')
    files_to_keep = []
    for i in range(100):
        files_to_keep.append(files[i*len(files)//100])
    frames = [Image.open(image) for image in files_to_keep]
    frame_one = frames[0]
    frame_one.save(f"{name}_prediction.gif", format="GIF", append_images=frames,
               save_all=True, duration=200, loop=0)
    

# visualize(states)