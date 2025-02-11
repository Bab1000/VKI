# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:10:27 2024

@author: mendez
"""


import os  # To create folders an delete them
import numpy as np

# Part 1: Get the Data
import urllib.request


FOLDER='/home/jpe/VKI/Courses/SP/Mendez_part/Tutorial_5_2D_Cylinder_Memory_Saving'
if not os.path.exists(FOLDER):
    os.mkdir(FOLDER)


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


# Plot a flow Field


k=2
FOLDER='/home/jpe/VKI/Courses/SP/Mendez_part/Tutorial_5_2D_Cylinder_Memory_Saving'+os.sep+'data' # Folder 
Name=FOLDER+os.sep+'Res%05d'%k+'.dat' # Check it out: print(Name)
# note: this is part of a modulo tutorial. To run the following 
# utility function, install modulo -> pip install modulo_vki
from modulo_vki.utils.others import Plot_Field_TEXT_Cylinder
Name_Mesh=FOLDER+os.sep+os.sep+'MESH.dat'
Name_FIG=FOLDER+os.sep+'Cylinder_Flow_snapshot_'+str(2)+'.png'

n_s,Xg,Yg,Vxg,Vyg,X_S,Y_S=Plot_Field_TEXT_Cylinder(Name,Name_Mesh,Name_FIG) 


# Construct the snapshot matrix
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
  # Get U component. Optional: reshape it as grid
  V_X=Dat[:,0]; #V_X_g=np.reshape(V_X,(n_y,n_x))
  # Get V component  Optional: reshape it as grid
  V_Y=Dat[:,1]; #V_Y_g=np.reshape(V_Y,(n_y,n_x))
  # Assign the vector like snapshots to snapshot matrices
  D_U[:,k]=V_X; D_V[:,k]=V_Y  
  print('Loading Step '+str(k+1)+'/'+str(n_t)) 

#% Store the snapshot matrices for the other exercises
np.savez('Snapshot_Matrices',D_U=D_U,D_V=D_V,Xg=Xg,Yg=Yg,t=t,X_S=X_S,Y_S=Y_S)


import matplotlib.pyplot as plt

# here is how you would make the plot of a certain snapshot:
fig, ax = plt.subplots(figsize=(5, 3)) # This creates the figure
k=11 # number of the snapshot to plot
U= D_U[:,k].reshape(n_y,n_x); V= D_V[:,k].reshape(n_y,n_x)     
Magn=np.sqrt(U**2+V**2)
CL=plt.contourf(Xg,Yg,Magn.T,levels=np.arange(0,18,2))
STEPx=1;  STEPy=1 # in case you want to jump some arrows
plt.quiver(Xg[::STEPx,::STEPy],Yg[::STEPx,::STEPy],\
               Vxg[::STEPx,::STEPy],Vyg[::STEPx,::STEPy],color='k')
plt.rc('text', usetex=True)      
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
fig.colorbar(CL,pad=0.05,fraction=0.025)
ax.set_aspect('equal') # Set equal aspect ratio
ax.set_xlabel('$x[mm]$',fontsize=13)
ax.set_ylabel('$y[mm]$',fontsize=13)
#ax.set_title('Tutorial 2: Cylinder Wake',fontsize=12)
ax.set_xticks(np.arange(0,70,10))
ax.set_yticks(np.arange(-10,11,10))
ax.set_xlim([0,50])
ax.set_ylim(-10,10)
circle = plt.Circle((0,0),2.5,fill=True,color='r',edgecolor='k',alpha=0.5)
plt.gcf().gca().add_artist(circle)
plt.tight_layout()
plt.show()


    
    
