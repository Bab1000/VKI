# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 12:05:41 2021

@author: mendez
"""

# In this exercise we use Finite Impulse Response filters
# to remove frequencies from a signal.

# The concepts to show are the following:
#
# 1. An ideal filter is not realizable. You shall not just zero frequencies
# 2. How to use Scipy to generate FIR filter's impulse response
# 3a. Show that the filter can be done via multiplication in freq domain 
#     or convolution in time domain. Show first this using python's functions.
# 3b. Show the same as above using linear algebra.
# 4. Show the filtfilt operation to compensate for delays.
# 5. Show how to use an IIR Filter

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import firwin


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


#%% Compute the FFT and and construct ideal Filter

# Define the cut off frequency (to play later)
f_c=50

# Step 1: Compute the  Fourier Transform of this signal
# We want to keep only frequencies below 50
# Discrete fast fourier transform
Signal_FFT = np.fft.fft(Signal)/np.sqrt(n_t) # Compute the DFT
Freqs=np.fft.fftfreq(len(Signal))*fs # Compute the frequency bins
plt.plot(Freqs,np.abs(Signal_FFT))
##plt.xlim([0,fs/2])

#%% Construct the Ideal Filter (usually not a good idea!)
Indices=np.argwhere(np.abs(Freqs) < f_c)
H_Transf=np.zeros(Freqs.shape)
H_Transf[Indices]=1;

# The impulse response of this signal:
h_impulse_ID=np.fft.fftshift(np.fft.ifft(H_Transf)*np.sqrt(n_t))
plt.plot(h_impulse_ID,'k--')    
    

# We could then do the filtering and invert
Signal_FFT_Filter=H_Transf*Signal_FFT
Signal_Filtered=np.fft.ifft(Signal_FFT_Filter)*np.sqrt(n_t)

#%% Construct the FIR Filter using Windowing Methods
N_O=int(len(Signal)/10) # this is the order of the filter
h_impulse = firwin(N_O, 50/fs*2, window='hamming')

# Determine the Transfer function of the filter.
# First we need to zero padd it!
n_zeros=len(Signal)-N_O
if n_zeros % 2 == 0:
    pass # Even number: we put equal zeros on right and left
    zero_left=int(n_zeros/2)
    zero_right=int(n_zeros/2)
else:
    pass # Odd number: we put one extra zero on the right
    zero_left=int((n_zeros-1)/2)
    zero_right=int(zero_left+1)

# for the  modulus we do not care about how the padding is done
h_impulse_padd=np.pad(h_impulse, (zero_left, zero_right), 'constant', constant_values=(0))

H_Transfer_actual=np.fft.fft(h_impulse_padd) # Compute the Transf Function DFT
# Compare the transfer function
H_Transfer_Abs=(np.abs(H_Transfer_actual))
plt.plot(Freqs[0:int(len(Freqs)/2)],H_Transfer_Abs[0:int(len(Freqs)/2)])
plt.xlim([0,90])

# There are now several methods to use the impulse response in the  frequency  domain
# See for example this https://scipy-cookbook.readthedocs.io/items/ApplyFIRFilter.html
# The crude way would be:
Sign_Filt_1=np.convolve(Signal,h_impulse,mode='same')
# Then, since this introduces a shift, we can filter it also backward
Sign_Filt_1b=np.flipud(np.convolve(np.flipud(Sign_Filt_1),h_impulse,mode='same'))

# Apply the filter from the diffence equation
Sign_FIR_DE=signal.lfilter(h_impulse, 1,Signal) # Use recursive formula 
plt.plot(t_k,Sign_Filt_1,'b-');plt.plot(t_k,Sign_FIR_DE,'r--')



# A Build in operation to cancel the phase shift of FIRs
from scipy.signal import filtfilt
Sign_Filt_2=filtfilt(h_impulse,1,Signal)


#%% Perform the Filtering using an IIR:
b, a = signal.butter(11, f_c, 'lp', analog=False,fs=fs)
Signal_IIR_Filtered=signal.lfilter(b,a,Signal) # Use recursive formula 

# Check the (digital) frequency response
#w_IIR, H_IIR = signal.freqz(b,a)
#plt.plot(w_IIR,abs(H_IIR),'k--')
# Check the impulse response (approximation)
# To approximate the impulse response we use this on an impulse
delta=np.zeros(t_k.shape); delta[0]=1 # Create Delta function
h_IIR=signal.lfilter(b,a,delta) # Use recursive formula on Delta function



# Compare the results (close to borders!)

fig, ax = plt.subplots(figsize=(9,5)) # Create Signal Noisy and Clean
plt.plot(t_k,Signal,label='Original')
plt.plot(t_k,Signal_Clean,label='Ideal Clean') # Show the clean (ideal signal)
#plt.plot(t_k,Sign_Filt_1,label='2 Convs')
#plt.plot(t_k,Sign_Filt_2,label='filtfilt')
plt.plot(t_k,Signal_Filtered,label='idealized')
plt.plot(t_k,Signal_IIR_Filtered,label='IIR')

plt.xlabel('$t_k[s] $',fontsize=12)
plt.ylabel('Signal',fontsize=12)
# plt.title('Eigen_Function_Sol_N',fontsize=18)
plt.xlim([0,0.2])
plt.ylim([-0.5,1]) 

plt.tight_layout()
plt.legend(fontsize=20,loc='upper left',ncol=3)
plt.show()
plt.savefig('Results_Filter_Border_1.pdf', dpi=100)      
plt.close(fig)




#%% We skip this part in class. Too much Algebra.
# If you are interested in knowing more, I am always ready to discuss Linear Algebgra

# What I am showing hereinafter is an interesting property, I used many times
# in my processing: I here show that the Fourier modes are eigenfunctions of the 
# Convolution operation and show the the associated eigenvalues are the entries of the 
# transfer function. This motivated me to develop a new decomposition.
# It might motivate you as well; if it doesn't, no worries. 
# You can still pass this exam.

# We now show that the transfer function is the set of eigenvalues of the 
# convolution matrix. Show that Psi_F H conj(Psi_F) is your Conv Matrix!
# First construct the Fourier Matrix in the fast way:
Psi_F=np.fft.fft(np.eye(n_t))
# diagonal matrix with the H
H_D=np.diag(H_Transfer_actual)
# If the calculation is correct, the convolution matrix should be:
C_conv=np.linalg.multi_dot([Psi_F,H_D,np.conj(Psi_F)])
# In each row of this, you should now see the impulse response shifted!
plt.plot(t_k,C_conv[1,:])
plt.plot(t_k,C_conv[10,:])
plt.plot(t_k,h_impulse_padd)










