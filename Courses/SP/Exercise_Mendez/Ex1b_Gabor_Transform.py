# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 15:13:28 2023

@author: admin
"""


# Exercise on Gabor Transform and Continuous Wavelet

import numpy as np
from numpy.fft import fft, fftshift
import matplotlib.pyplot as plt  # This is to plot things


plt.rc('text', usetex=True)  # This is Miguel's customization
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)


T = 1 # Duration of the signal
f_s=1000 # Sampling frequency in Hz
dt=1/f_s # Sampling Frequency
ta=np.arange(0,T,dt) # time axis 
n_t=len(ta) # Number of points in time
tb=np.arange(0,T/4,dt) # time axis 
  
frequencies = [4, 30, 60, 90]
y1a, y1b = np.sin(2*np.pi*frequencies[0]*ta), np.sin(2*np.pi*frequencies[0]*tb)
y2a, y2b = np.sin(2*np.pi*frequencies[1]*ta), np.sin(2*np.pi*frequencies[1]*tb)
y3a, y3b = np.sin(2*np.pi*frequencies[2]*ta), np.sin(2*np.pi*frequencies[2]*tb)
y4a, y4b = np.sin(2*np.pi*frequencies[3]*ta), np.sin(2*np.pi*frequencies[3]*tb)
 
# Assembly the two signals in the time domain
s1 = y1a + y2a + y3a + y4a
s2 = np.concatenate([y1b, y2b, y3b, y4b])
 

# Compute the spectra of the two
n_omega=len(s1) # Number of frequency bins
df=f_s/n_omega
Freq_d=np.linspace(-f_s/2,f_s/2-df,n_omega)
Freq_d2=fftshift(np.fft.fftfreq(n_omega))*f_s # This is np tool for freq axis


# We now proceed with the Gabor Transform of the Signal 2
# First we define the width of our Gabor Kernel, 
# and plot it in the middle of the domain
ALPHAS=np.array([5,500,5000])

for i in range(0,len(ALPHAS)):
 alpha=ALPHAS[i]
 time=ta-0.5;
 G_ex=(2*alpha*np.pi)**0.25*np.exp(-np.pi*alpha*time**2)
 LAB_string='$\\alpha='+str(alpha)+'$'
 plt.plot(ta,G_ex,label=LAB_string)

plt.legend()
plt.show()


# Before goind for the set of convolutions, we define a time shift discretizations
Taus=np.linspace(0,1,200) # This defines the time resolution of your grid
################### CHANGE THIS ###########################################
alpha=500 # try 5 500 or 5000
##########################################################################

# We will also create a temporary folder to store the 
G_Grid=np.zeros([n_t,len(Taus)])

import os

GIFNAME = 'Gabor_alpha_'+str(alpha)+'.gif'
Fol_Out = 'Gif_Images_temporary'
if not os.path.exists(Fol_Out):
    os.mkdir(Fol_Out)
    
# We export one image every NS
NS=10
INDEX_E=np.arange(1,len(Taus),NS)
counter=0# We initialize a value for the saving.

for k in range(1,len(Taus)):
  # Construct the Gabor Kernel at each shift
  G_Ker=(2*np.pi*alpha)**0.25*np.exp(-np.pi*alpha*(ta-Taus[k])**2)
  # Multiply Kernel and Signal
  F_G=G_Ker*s2
  # Compute the Fourier Transform of this
  s_g_hat=fftshift(fft(F_G)) # Fourier Transform 1
  # Store the result in the G_Transform
  G_Grid[:,k]=np.abs(s_g_hat)
  # Export every NS images: if the k is in INDEX_E, print!
  if (k in INDEX_E):
      fig, axarr = plt.subplots(figsize=(6,8))
      plt.subplot(3,1,1)
      plt.plot(ta,s2)
      plt.plot(ta,G_Ker/np.max(G_Ker),color='r',linewidth=2.0)
      plt.xlabel('t [s]',fontsize=18)
      plt.title('Signal and Kernel',fontsize=18)
      plt.subplot(3,1,2)
      plt.plot(ta,F_G)
      plt.xlabel('t [s]',fontsize=22)
      plt.title('Windowed Signal',fontsize=18)
      plt.subplot(3,1,3)
      plt.plot(Freq_d,G_Grid[:,k])
      plt.xlabel('f [Hz]',fontsize=22)
      plt.xlim([-160,160])
      plt.title('Spectra of Windowed Signal',fontsize=18)
      plt.tight_layout()
      NameOUT=Fol_Out + os.sep + 'Im%03d' % (counter) + '.png'
      plt.savefig(NameOUT, dpi=100)
      counter=counter+1
      plt.close(fig)
  print('Image n ' + str(k) + ' of ' + str(n_t))

plt.close()


# Assembly the GIF
import imageio  # This used for the animation

images = []
LENGHT = len(INDEX_E)
counter=0

for k in range(1, LENGHT, 1):
    MEX = 'Preparing Im ' + str(k) + ' of ' + str(LENGHT)
    print(MEX)    
    FIG_NAME = Fol_Out + os.sep + 'Im%03d' % (counter) + '.png'
    images.append(imageio.imread(FIG_NAME))
    counter=counter+1

# Now we can assembly the video and clean the folder of png's (optional)
imageio.mimsave(GIFNAME, images, duration=0.2)
import shutil  # nice and powerfull tool to delete a folder and its content

shutil.rmtree(Fol_Out)

# We now plot the time-frequency plane of the Gabor Transform

fig = plt.subplots(figsize=(6,4))
plt.contourf(Taus,Freq_d,G_Grid)
plt.title('Gabor Transform',fontsize=18)
plt.ylim([-160,160])
plt.ylabel('f [Hz]',fontsize=18)
plt.xlabel('t [s]',fontsize=18)
plt.tight_layout(pad=1, w_pad=0.5, h_pad=1.0)
NameOUT='Gabor_T_alpha_'+str(alpha)+'.png'
plt.savefig(NameOUT, dpi=100)






