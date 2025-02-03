# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:10:27 2024

@author: mendez
"""


#%% Part 1: Get the Data
import urllib.request
print('Downloading PIV data TR PIV Cyl...')
url = 'https://osf.io/47ftd/download'
urllib.request.urlretrieve(url, 'Ex_5_TR_PIV_Cylinder.zip')
print('Download Completed! I prepare data Folder')
# Unzip the file 
from zipfile import ZipFile; import shutil
String='Ex_5_TR_PIV_Cylinder.zip'
zf = ZipFile(String,'r'); 
zf.extractall('./DATA_CYLINDER'); zf.close()
shutil.move('DATA_CYLINDER', FOLDER+os.sep+'data') # rename the data flolder to FOLDER
os.remove(String) # Delete the zip file with the data 
print('Data set unzipped and ready ! ')

n_t=13200; Fs=3000; dt=1/Fs 
t=np.linspace(0,dt*(n_t-1),n_t) # prepare the time axis# 


#%% Plot a flow Field
import os
import numpy as np

FOLDER='Tutorial_5_2D_Cylinder_Memory_Saving'+os.sep+'data' # Folder 
Name=FOLDER+os.sep+'Res%05d'%k+'.dat' # Check it out: print(Name)
Name_Mesh=FOLDER+os.sep+'MESH.dat'


#%% Construct the snapshot matrix
# number of points in x and y
n_x,n_y=np.shape(Xg)
# Prepare the time axis
n_t=13200; Fs=3000; dt=1/Fs 
t=np.linspace(0,dt*(n_t-1),n_t) # prepare the time axis# 
# number of velocity points 
#(from Plot_Field_TEXT_Cylinder)
nxny=int(n_s/2) 

# Prepare the snapshot matrices U and V
D_U=np.zeros((n_s//2,n_t),dtype="float32") 
D_V=np.zeros((n_s//2,n_t),dtype="float32") 

for k in range(0,n_t):
  # Name of the file to read
  Name=FOLDER+os.sep+'Res%05d'%(k+1)+'.dat' 
  # Read data from a file
  # Here we have the two colums
  DATA = np.genfromtxt(Name,usecols=np.arange(0,2),max_rows=nxny+1) 
  Dat=DATA[1:,:] # Remove the first raw with the header
  # Get U component and reshape as grid
  V_X=Dat[:,0]; V_X_g=np.reshape(V_X,(n_y,n_x))
  # Get V component and reshape as grid
  V_Y=Dat[:,1]; V_Y_g=np.reshape(V_Y,(n_y,n_x))
  # Assign the vector like snapshots to snapshot matrices
  D_U[:,k]=V_X; D_V[:,k]=V_Y  
  print('Loading Step '+str(k+1)+'/'+str(n_t)) 

#% Store the snapshot matrices for the other exercises
np.savez('Snapshot_Matrices',D_U=D_U,D_V=D_V,Xg=Xg,Yg=Yg,t=t,X_S=X_S,Y_S=Y_S)