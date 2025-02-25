# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 18:25:26 2023

@author: mendez
"""

# In this exercise we use Finite Impulse Response filters
# to remove frequencies from a signal.

# The concepts to show are the following:
#
# 1. Understand recursive relations in IIR filters
# 2. Use the signal processing toolbox to generate the filter coefficients
# 3  Implement the recursive relation 'by hand'
# 4. Implement the recusrive relation with python's function lfilter

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


plt.rc('text', usetex=True)      # This is for plot customization
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


#%% Construct Signal
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


fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
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


#%% Step 2: Check out the spectral of this signal

# Finally, we do it with python's fft (much faster)
u_hat_fft=(np.fft.fft(Signal)) # this is not the nice scaling.. left as exercise :)
freqs_fft=freq = np.fft.fftfreq(n_t, d=dt) # Check these out !


# plot the results 
plt.figure()
plt.plot(freqs_fft,np.abs(u_hat_fft))
plt.xlabel('freqs $[Hz]$',fontsize=24)
plt.ylabel('u hat Abs',fontsize=24)
plt.tight_layout()
plt.savefig('Absolute Values_shifted_and_f.png', dpi=100)      
plt.show()

# a more canonical way to plot this:
plt.figure()
plt.plot(freqs_fft[:n_t//2],2*np.abs(u_hat_fft[:n_t//2])/n_t)
plt.xlabel('freqs $[Hz]$',fontsize=24)
plt.ylabel('u hat Abs',fontsize=24)
plt.xscale('log',base=10);plt.yscale('log',base=10)
plt.tight_layout()
plt.savefig('Positive_Part_only.png', dpi=100)      
plt.show()

    
#%% Step 3: Compute the filter coefficient for a low pass with f_c=50Hz
# Here we use the butter implementation from scipy
from scipy import signal
f_c=200 # cut off frequency in Hz
b, a = signal.butter(3, f_c, 'lp', analog=False,fs=fs) # Filter coeffs
# This creates the transf function obj
sys = signal.TransferFunction(b, a, dt=1) 
w, mag, phase = sys.bode() # Note that the result is in db by default !

# Check out the digital transfer functions.
# The digital frequencies are (take -pi to pi):
dtheta=2*np.pi/n_t
theta_n=np.linspace(0,n_t-1,n_t)*dtheta-np.pi


# Order 3
b3, a3 = signal.butter(3, f_c, 'lp', analog=False,fs=fs) # Filter coeffs
sys = signal.TransferFunction(b3, a3, dt=1) 
w3, mag3, phase3 = sys.bode(w=theta_n) # Note that the result is in db by default !


# Order 7
b7, a7 = signal.butter(7, f_c, 'lp', analog=False,fs=fs) # Filter coeffs
sys = signal.TransferFunction(b7, a7, dt=1) 
w7, mag7, phase7 = sys.bode(w=theta_n) # Note that the result is in db by default !

# Order 11
b11, a11 = signal.butter(11, f_c, 'lp', analog=False,fs=fs) # Filter coeffs
sys = signal.TransferFunction(b11, a11, dt=1) 
w11, mag11, phase11 = sys.bode(w=theta_n) # Note that the result is in db by default !


# Let's plot the magnitudes and the phases (in a linear plot!)
fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
plt.plot(w3,10**(mag3/10),label='O3')
plt.plot(w7,10**(mag7/10),label='O7')
plt.plot(w11,10**(mag11/10),label='O11')
plt.xlim([-np.pi,np.pi])
plt.legend(fontsize=20)
plt.xlabel('digital freqs $[rad/s]$',fontsize=24)
plt.ylabel('H Abs',fontsize=24)
plt.tight_layout()
plt.savefig('H_IIR.png', dpi=100)      
plt.show()

# The plot in the classic log vs decibel form are the following:
fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
plt.semilogx(w3/np.pi*fs/2,mag3,label='O3')
plt.semilogx(w3/np.pi*fs/2,mag7,label='O7')
plt.semilogx(w3/np.pi*fs/2,mag11,label='O11')
plt.legend(fontsize=20)
plt.xlabel('f $[Hz]$',fontsize=24)
plt.ylabel('H Abs',fontsize=24)
plt.tight_layout()
plt.savefig('H_IRR_log.png', dpi=100)      
plt.show()
    
    
# The plot in the classic log vs decibel form are the following:
fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
plt.semilogx(w3/np.pi*fs/2,phase3,label='O3')
plt.semilogx(w3/np.pi*fs/2,phase7,label='O7')
plt.semilogx(w3/np.pi*fs/2,phase11,label='O11')
plt.legend(fontsize=20)
plt.xlabel('f $[Hz]$',fontsize=24)
plt.ylabel('H phase',fontsize=24)
plt.tight_layout()
plt.savefig('H_IRR_log_ph.png', dpi=100)      
plt.show()


#%% Apply the filter with a recureive formula (only Filter 33)


# Let us consider the first filter.
y=np.zeros_like(Signal)
u=Signal # This is the input signal from the slide's notation

# This should be the general recursrive form: 
      
for k in range(len(y)):
   if k==0:
     print('k='+str(k))
     y[k]=b3[0]*u[k]  # First Step
   if k==1:
     print('k='+str(k))
     y[k]=b3[0]*u[k]+b3[1]*u[k-1] -a3[1]*y[k-1]  # Second Step
   elif k==2:
     print('k='+str(k))
     y[k]=b3[0]*u[k]+b3[1]*u[k-1]+b3[2]*u[k-2]-a3[1]*y[k-1]-a3[2]*y[k-2]
   else:
     print('k='+str(k))
     y[k]=b3[0]*u[k]\
           +b3[1]*u[k-1]\
           +b3[2]*u[k-2]\
           +b3[3]*u[k-3]\
           -a3[1]*y[k-1]\
           -a3[2]*y[k-2]\
           -a3[3]*y[k-3]
     
     
# Using Python's lfilter:
y_p=signal.lfilter(b3,a3,u) # Use recursive formula 
    

# Filter order 3: Python vs manual
fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
plt.plot(t_k,y,label='Manual')
plt.plot(t_k,y_p,label='Scipy')
plt.legend(fontsize=20)
plt.xlabel('$t_k [s]$',fontsize=24)
plt.ylabel('$y$',fontsize=24)
plt.tight_layout()
plt.savefig('u_filters_h_vs_p.png', dpi=100)      
plt.show()
    
# We can compare the spectra filtered vs nonfiltered:
u_hat_fft_filt=(np.fft.fft(y_p)) # spectra of the filtered signal

    
# a more canonical way to plot this:
plt.figure()
plt.plot(freqs_fft[:n_t//2],2*np.abs(u_hat_fft[:n_t//2])/n_t,label='Original')
plt.plot(freqs_fft[:n_t//2],2*np.abs(u_hat_fft_filt[:n_t//2])/n_t,label='Filtered')
plt.legend(fontsize=20)
plt.xlabel('freqs $[Hz]$',fontsize=24)
plt.ylabel('u hat Abs',fontsize=24)
plt.xscale('log',base=10);plt.yscale('log',base=10)
plt.tight_layout()
plt.savefig('Filtered_signal.png', dpi=100)      
plt.show()
    
    
    
    
    
    
    
    



