# -*- coding: utf-8 -*-
"""
Created on Fri Nov  24 08:55:41 2023

@author: mendez
"""

# In this exercise we play with the Fourier Transform

# The concepts to show are the following:
#
# 1. The Fourier Transform is a set of dot products
# 2. The DFT assumes that your signal is Periodic! 
# 3. The Fourier Matrix is the most beautiful matrix ever.

# 4. Bonus: The Fourier has no information about time ! When things happen?
# this is a separate file.


import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)      # This is for plot customization
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


#%% Construct Signal

# Note: The normalization we are following here is to have 
# the Fourier Matrix of unitary length. This is not the normalization you
# would use to preserve the amplitude of the sinusoids.

tn = 0.58 # [s] end time
fs = 1000 #[Hz] sampling frequency
Nt = int(fs*tn)

fext1 = 20 # [Hz] signal frequency
fext2 = 5 # [Hz] signal frequency


Noise_scale=0.1 # pick a value of the order of 0.3 ?


#time signal
t = np.linspace(0, tn, num=Nt) 
#sine signal 
y = 2*np.sin(2*np.pi*fext1*t) + 3*np.sin(2*np.pi*fext2*t)+\
    np.random.normal(scale=Noise_scale, size=Nt)

#plot signal 
plt.figure()
plt.plot(t,y)
plt.xlabel('time [s]')
plt.ylabel('Signal amplitude')
plt.show()

print('The length of the initial time signal y is %i' %(Nt))



#%% DFT version 1 : many dot products ! 

Nf = Nt # frequency sampling try [<Nt =Nt >Nt]

# radians vector
omega_DTFT = np.linspace(-np.pi,np.pi,num = Nf, endpoint = False)
#tranformation to frequency
f = omega_DTFT*fs/(2*np.pi)

# initialize DTFT vector
y_hat_DFTF = np.zeros(len(f), dtype=complex)

# important observation: np.dot has no conjugation!!

#compute the DFTF as a dot product between the signal and the fourier basis
for jj in range(len(f)):
    y_hat_DFTF[jj] = np.dot(y,np.exp(-1j*omega_DTFT[jj]*np.arange(Nt)))/np.sqrt(Nt)
           
plt.figure()
plt.plot(f,np.abs(y_hat_DFTF))
plt.xlabel('frequency [Hz]')
plt.ylabel('DFTF Abs')
plt.show()

plt.figure()
plt.plot(f,np.angle(y_hat_DFTF))
plt.xlabel('frequency [Hz]')
plt.ylabel('DFTF Angle')
plt.show()

#%% DFT version 2 : The Fourier matrix (for the moment only if Nt=Nf )
if Nt == Nf:
    
    #intitialize
    PSI_F = np.zeros((len(y), len(omega_DTFT)), dtype=complex)
    w=np.exp(1j*2*np.pi/Nt)
    #build the Fourier matrix
    for ii in range(len(y)):   
        for jj in range(len(omega_DTFT)):
           # PSI_F[ii,jj] =  np.exp(1j*omega_DTFT[jj]*ii)/np.sqrt(Nt)
            PSI_F[ii,jj] =  w**(ii*jj)/np.sqrt(Nt)
    
    #compute complex conjugate
    CJ= np.conj(PSI_F)    
    
    #multiplication with fourier matrix
    y_hat_DFT = np.dot(CJ,y)
            
    PSI_F2=np.conj(np.fft.fft(np.eye(len(y))))/np.sqrt(len(y))
    # These are some nice check: How to show that the Fourier basis is orthonormal?
    #d_hat=np.dot(np.conj(PSI_F),y)
    #d=np.dot(PSI_F,d_hat)
     
    # Identity=(np.abs(np.dot(np.conj(PSI_F),PSI_F)))
    # Identity[0,0]
    # Identity=(np.abs(np.dot(np.conj(PSI_F2),PSI_F2)))
    # Identity[0,0]
    
    #dispaly PsiF matrix
    fig, ax = plt.subplots()
    cax = ax.matshow(np.real(PSI_F))
    cbar = fig.colorbar(cax)
    ax.set_title('Real part of Fourier matrix')
    fig.show()
    
    fig, ax = plt.subplots()
    cax = ax.matshow(np.abs(np.dot(CJ.T,PSI_F)))
    cbar = fig.colorbar(cax)
    ax.set_title('Check matrix orthogonality')
    plt.show()

#%% DFT version 3: the FFT

# Discrete fast fourier transform
    DFFT = np.fft.fft(y)/np.sqrt(Nt)
    
    plt.figure()
    plt.plot(f,np.angle(y_hat_DFTF),color='b',label='DFTF')
    plt.plot(f,np.angle(np.fft.fftshift(DFFT)),color='r', linestyle='', marker='v',markersize=5,label='DFFT')
    plt.plot(f,np.angle(np.fft.fftshift(y_hat_DFT)),color='g', linestyle='--',label='DFT')
    plt.xlabel('frequency [Hz]')
    plt.legend()
    plt.xlim([-19,19])
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.plot(f,np.abs(y_hat_DFTF),color='b',label='DFTF')
    plt.plot(f,np.abs(np.fft.fftshift(DFFT)),color='r', linestyle='', marker='v',markersize=5,label='DFFT')
    plt.plot(f,np.abs(np.fft.fftshift(y_hat_DFT)),color='g', linestyle='--',label='DFT')
    plt.xlabel('frequency [Hz]')
    plt.legend()
    plt.xlim([-19,19])
    plt.grid()
    plt.show()



