# -*- coding: utf-8 -*-
"""
Created on Jan 23 15:03:41 2024

@author: mendez
"""

import numpy as np
import matplotlib.pyplot as plt


#%% Preamble: customization of matplotlib
# Configuration for plots
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


#%% Harmonic function definition
def harmo(X):  # The Function to Minimize  
    """Function to be minimized"""  # new        
    C=X[0]*np.sin(4*X[0])+1.1*X[1]*np.sin(2*X[1])
    return C 

#%% Function implementing the Simulated Annealing
def Simulated_A(func,X0,X_Bounds,n_iter,R,t_I,a=1):
  # Uses simulated annealing to minimize func, starting from X0
  # and running for n_iterations within the domain X_Bounds and using
  # a box +-R for the new search. The number of iterations is n_iter.
  # The initial temperature is t_I and a is the exponent in the cooling
  # law. The output is the full story of X_k, func(X_k) and the best result.
  n_dim=np.shape(X_Bounds)[0] 
  # Initialize the next vector  
  X_next=np.zeros(n_dim)
  # Initialize the new vector
  X_best=X0;
  # Prepare the story of attempt and evaluations
  X_k,Cost_k =np.zeros((n_dim,n_iter+1)),np.zeros(n_iter+1)
  
  # Run the first evaluation
  X_k[:,0]=X_best
  Cost_k[0]=func(X_best)
  F_best=Cost_k[0] # Define the initial best
  print('Current Best is X='+str(X_best))
  
  # Loop over the iterations  
  for k in range(1,n_iter):    
    # Pick a random vector in the given range
    for j in range(n_dim):
       X_next[j]=np.random.uniform(low=np.max([X_k[j,k-1]-R,X_Bounds[j][0]]),
                                   high=np.min([X_k[j,k-1]+R,X_Bounds[j][1]]))        
    # Evaluate function in the new point
    F_next=func(X_next)
    # Check if this solution is better (if so, store it)
    if F_next<F_best:
       X_best=X_next
       F_best=func(X_best)
       print('Change the current Best to X='+str(X_best))
    #------------------------- Annealing Process-------------------------   
    # Get 'Temperature'
    T=t_I/(k**a+1) 
    # Get the difference between current and best (could be 0!)
    diff=np.abs(F_next-F_best)/F_best
    # Compute the metropolis acceptance
    Metropolis=np.exp(-diff/T)
    Rand_Change=np.random.rand()<Metropolis
    # Tell us if the change occurred randomly
    if Rand_Change:
       print('Random Change Possible! High Temperature!')  
    # We might do a replacement even if it's 'Bad'
    if F_next<Cost_k[k-1] or Rand_Change: 
     # Go for the new option   
     X_k[:,k]=X_next
     Cost_k[k]=F_next
    else:
     # Stay where you are   
     X_k[:,k]=X_k[:,k-1]
     Cost_k[k]=Cost_k[k-1]    
    
  return X_k,Cost_k,X_best    
      






#%% Make animation of the RS History
import os
from matplotlib import cm

import imageio


def Anim(X_k,Cost_k,func,x_1m,x_1M,x_2m,x_2M,n_p,Name_Video):
   # This function makes an animation of the single search 
   # The X_k is the story of the optimizer.
   # func is the function optimized
   # x_1m,x_1M,x_2m,x_2M,n_p is the same input for the plot func
   # Name_Video is the name of the gif that will be exported
    
   # Temporary Folder
   FOLDER='Temp'
   if not os.path.exists(FOLDER):
    os.makedirs(FOLDER) 
   
   #%% Prepare the Contour
   # Create the vectors for the grid  
   x = np.linspace(x_1m, x_1M, n_p)
   y = np.linspace(x_2m, x_2M, n_p)
   X, Y = np.meshgrid(x, y) # Create grid
   COST=np.zeros((n_p,n_p)) # Initialize the cost func
   # Evaluate the cost function
   for i in range(0,len(x)):
     for j in range(0,len(x)):
      XX=np.array([X[i,j],Y[i,j]]) # Get Point Loc   
      COST[i,j]=func(XX) # Interrogate the func   
   
   # Find the approximate location of the minima
   obb=COST.min()
   ID=np.where(COST == obb) 
   
   # Get number of iterations
   n_iter=np.shape(X_k)[1]
   plt.ioff()

   for k in range(n_iter):
    print('Iteration '+str(k)+'/'+str(n_iter))   
    fig= plt.figure(figsize=(10, 4)) # This creates the figure
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    ax1.contourf(X, Y, COST,cmap=cm.coolwarm, alpha=0.8) 
    ax1.plot(X[ID],Y[ID],'wo',markersize=5)
    ax1.plot(X_k[0,0:k],X_k[1,0:k],'ko:')
    
    ax1.set_xlim([x_1m, x_1M])
    ax1.set_ylim([x_2m, x_2M])           
  
    plt.xticks(np.arange(x_1m,x_1M,1))
    plt.yticks(np.arange(x_2m,x_2M,1)) 
    
    ax2.plot(np.linspace(0,k,k),Cost_k[0:k],'ro:')
    ax2.set_ylim([1.1*np.min(Cost_k),1.2*np.max(Cost_k)])
    ax2.set_yticks([])
                   
                   
    #ax2.set_xlim([0,k])
    ax2.set_xticks([])
    
    
    plt.title("Iteration "+ str(k))
    plt.savefig(FOLDER+'/'+'Step'+str(k)+'.png')
    plt.close('all')
  
   # Make a Gif 1
   GIFNAME=Name_Video
   images=[]    
   for k in range(n_iter):
     MEX= 'Preparing Im '+ str(k)+' of ' + str(n_iter-1)
     print(MEX)
     FIG_NAME=FOLDER+'/'+'Step'+str(k)+'.png'
     images.append(imageio.imread(FIG_NAME))
   
   imageio.mimsave(GIFNAME, images,duration=0.5)
   import shutil  # nice and powerfull tool to delete a folder and its content
   shutil.rmtree(FOLDER)


#%% Compact Code
X_Bounds=[(-8,8),(-8,8)]  # Boundaries for x1 and x2
n_iter=100
X0=np.array([-2,0])
t_I=harmo(X0) # Initial temperature
R=0.2 # Width of the exploration

X_k,Cost_k,X_best =Simulated_A(harmo,X0,X_Bounds,n_iter,R,t_I,1)

x_1m,x_1M,x_2m,x_2M,n_p=-8,8,-8,8,200
Anim(X_k,Cost_k,harmo,x_1m,x_1M,x_2m,x_2M,n_p,'harmo_SA.gif')




