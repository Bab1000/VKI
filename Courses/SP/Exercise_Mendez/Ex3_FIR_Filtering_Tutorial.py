# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 09:29:30 2023

@author: admin
"""

# This tutorial is about FIR filters.
# We create the same signal as last time, then we :
    
# 1. create the impulse response of FIR filters and study the bode plots
# 2. Apply the FIR using recursive relations, convolutions and ffts. 
# 3  Implement them in a non causal way (back-forth)
# 4. We close with a bonus by Miguel: the MRA analysis of a signal. The proto-wavelet!


import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import firwin


plt.rc('text', usetex=True)      # This is for plot customization
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


#%% Step 0: generate the same signal as last time
# this is a copy-paste from the previous exercise
# Fist we construct a signal that has several harmonics plus some noise
n_t=np.power(2,12) # number of time steps
fs=3000; dt=1/fs#  Sampling frequency in Hz
# Time Axis
t_k=np.arange(0,n_t*dt,dt)#Time Discretization

# We construct the signal as the sum of various contributions
F_1=40; # frequency of the signal
sigma_1=0.2 # Standard deviation of the modulation
X_S_1=0.7; # Location of the modulation
Signal_1=np.sin(2*np.pi*F_1*t_k)*(np.exp(-(t_k-X_S_1)**2/(2*sigma_1**2)))
plt.plot(t_k,Signal_1)

# We construct the signal as the sum of various contributions
F_2=400; # frequency of the signal
sigma_2=0.1 # Standard deviation of the modulation
X_S_2=0.8; # Location of the modulation
Signal_2=np.sin(2*np.pi*F_2*t_k)*(np.exp(-(t_k-X_S_2)**2/(2*sigma_2**2)))
plt.plot(t_k,Signal_2)


# Add a very large scale 
F_3=0.5; # frequency of the signal
sigma_3=1 # Standard deviation of the modulation
X_S_3=0.7; # Location of the modulation
Signal_3=np.sin(2*np.pi*F_3*t_k)*(np.exp(-(t_k-X_S_3)**2/(2*sigma_3**2)))
plt.plot(t_k,Signal_3)

# Add also some random noise 
Noise=0.1*np.random.randn(len(Signal_2))

# Add all the three signals
plt.close()
Signal=Noise+Signal_1+Signal_2+Signal_3+0.2

# # We want to keep only frequencies below 50
# Ideally the clean signal should be:
Signal_Clean=Signal_1+Signal_3+0.2


fig, ax = plt.subplots(figsize=(9,5)) 
plt.plot(t_k,Signal,label='Original')
plt.plot(t_k,Signal_Clean,label='Ideal Clean') # Show the clean (ideal signal)

plt.xlabel('$t_k[s] $',fontsize=12)
plt.ylabel('Signal',fontsize=24)
# plt.title('Eigen_Function_Sol_N',fontsize=18)
plt.xlim([0,np.max(t_k)])
#plt.ylim([0,1.1]) 
plt.tight_layout()
plt.legend(fontsize=24)
plt.savefig('Signal_to_Filter.pdf', dpi=100)      
plt.close(fig)



#%% Step 1: create impulse response and bode plots

# Define the cut off frequency
f_c=350
# Define the order of the filter:
N_O=150 # a good rule of thumb: int(len(Signal)/10)
h_impulse = firwin(N_O, f_c/fs*2, window='hamming')

# Have a look at it:
plt.figure(11)
plt.plot(h_impulse,'k-')  
plt.xlabel('$k $',fontsize=12)
plt.ylabel('h[k]',fontsize=12)

# Check the transfer function of this  (compare this for the IIR tutorial)
dtheta=2*np.pi/n_t
theta_n=np.linspace(0,n_t-1,n_t)*dtheta-np.pi
sys = signal.TransferFunction(h_impulse,1,dt=1) 
w, mag, phase = sys.bode(w=theta_n) # Note that the result is in db by default !

    
# Let's plot the magnitudes and the phases (in a linear plot!)
fig, ax = plt.subplots(figsize=(9,5))
plt.plot(w,10**(mag/10),label=['order_'+str(N_O)]) 
plt.xlim([-np.pi,np.pi])
plt.legend(fontsize=20)
plt.xlabel('digital freqs $[rad/s]$',fontsize=24)
plt.ylabel('H Abs',fontsize=24)
plt.tight_layout()
plt.savefig(['H_FIR_NO_'+str(N_O)+ '.png'], dpi=100)      
plt.show()  


# The plot in the classic log vs decibel form are the following:
fig, ax = plt.subplots(figsize=(9,5)) 
plt.semilogx(w/np.pi*fs/2,mag,label=['order_'+str(N_O)])
plt.legend(fontsize=20)
plt.xlabel('f $[Hz]$',fontsize=24)
plt.ylabel('H Abs',fontsize=24)
plt.ylim([-150,1])
plt.tight_layout()
plt.savefig(['H_FIR_Log_NO_'+str(N_O)+ '.png'], dpi=100)      
plt.show()


#%% Step 2: Implement the FIR filter in many ways ! 

# There are many ways of applying a filter. Here I show you 4.
# For more methods have a look at https://scipy-cookbook.readthedocs.io/items/ApplyFIRFilter.html


u=Signal # just to keep the same notation as in the slides
# Way 1: recursive relation (causal)
y_rec=signal.lfilter(h_impulse,1,u) # Use recursive formula 

# Way 2: standard (causal) convolution
y_conv_nc=np.convolve(u,h_impulse,mode='same')


# Way 3: non-causal convolution
h_impulse_pad=np.pad(h_impulse, (N_O, 0), 'constant', constant_values=(0))    
y_conv=np.convolve(u,h_impulse_pad,mode='same')

# Way 4:in the frequency domain (non-causal and implying periodicity) 
y_fft=np.real(np.fft.ifft(np.fft.fft(h_impulse,len(u))*np.fft.fft(u)))



fig, ax = plt.subplots(figsize=(9,5)) 
plt.plot(u,label='O')
plt.plot(y_rec,'rs:',label='rec')
plt.plot(y_conv_nc,'ko:',label='conv nc')
plt.plot(y_conv,label='conv')
plt.plot(y_fft,label='fft')
plt.legend(fontsize=20)
plt.xlabel('k ',fontsize=24)
plt.ylabel('y',fontsize=24)
plt.tight_layout()
plt.savefig('FIR_Filtered.png', dpi=100)      
plt.show()



#%% Step 3: Implement the FIR filter in many ways ! 

# Manual approach to face cancellation:
# Forward filter
Sign_Filt_1=signal.lfilter(h_impulse,1,u)
# Backward filter
Sign_Filt_1b=np.flipud(signal.lfilter(h_impulse,1,np.flipud(Sign_Filt_1)))

from scipy.signal import filtfilt
Sign_Filt_2=filtfilt(h_impulse,1,Signal)


fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
plt.plot(u,label='O')
plt.plot(Sign_Filt_1b,'rs:',label='Manual')
plt.plot(Sign_Filt_2,'ko:',label='filtfilt')
plt.legend(fontsize=20)
plt.xlabel('k ',fontsize=24)
plt.ylabel('y',fontsize=24)
plt.tight_layout()
plt.savefig('manual_vs_filt_filt.png', dpi=100)      
plt.show()



# BONUS QUESTION: how to help the boundaries: what happens if we try to tape
# the signal to zero ? for example we could subtrack a line connecting the first 
# and the last point and then add it back ?

# Define line through the points 1 and end:
x_f=np.array([0,len(u)]); y_f=np.array([0.17,-0.6])
p_line=np.polyfit(x_f,y_f,1)
u_line=np.polyval(p_line,np.arange(0,len(u),1))
# The shift is :
u_shift=u-u_line


fig, ax = plt.subplots(figsize=(9,5)) 
plt.plot(u,label='Original')
plt.plot(u_line,'r',label='Line Shift')    
plt.plot(u_shift,'k',label='shifted')
plt.legend(fontsize=20)
plt.xlabel('k ',fontsize=24)
plt.ylabel('y',fontsize=24)
plt.tight_layout()
plt.savefig('Shift_line.png', dpi=100)      
plt.show()


# Now go for the filtering:
Sign_Filt_2_shift=filtfilt(h_impulse,1,u_shift)+ u_line 
    

fig, ax = plt.subplots(figsize=(9,5)) 
plt.plot(u,label='O')
plt.plot(Sign_Filt_2_shift,'rs:',label='Shifted Filt')
plt.plot(Sign_Filt_2,'ko:',label='filtfilt')
plt.legend(fontsize=20)
plt.xlabel('k ',fontsize=24)
plt.ylabel('y',fontsize=24)
plt.tight_layout()
plt.savefig('manual_vs_filt_filt.png', dpi=100)      
plt.show()



#%% Step 5: (Bounus)

# Here we perform Multi resolution analysis
from MRA_Functions_Mendez import MRA_SISO_FILT_FILT

F_V=np.array([10,100,600])
U_MRA,H_MRA=MRA_SISO_FILT_FILT(u,F_V,fs,200)

# Show the four scales:
fig, ax = plt.subplots(figsize=(9,5)) 
plt.plot(U_MRA[:,0],label='scale1') 
plt.plot(U_MRA[:,1],label='scale2')  
plt.plot(U_MRA[:,2],label='scale3')
plt.plot(U_MRA[:,3],label='scale4')  
plt.legend(fontsize=20)
plt.xlabel('k ',fontsize=24)
plt.ylabel('y',fontsize=24)
plt.tight_layout()
plt.savefig('identified_scales_MRA.png', dpi=100)      
plt.show()





