# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 17:17:27 2020

@author: mendez
"""

import numpy as np
import matplotlib.pyplot as plt



plt.rc('text', usetex=True)      # This is for plot customization
plt.rc('font', family='serif')
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

n_t=4096 # number of points
# Time discretization
k=np.arange(0,n_t,1);f_s=1000; dt=1/f_s; t_k=k*dt
# Prepare three terms
u_1=2*np.sin(2*np.pi*20*t_k)*np.exp(-(t_k-1)**2/0.05)
u_2=1*np.sin(2*np.pi*90*t_k)*np.exp(-(t_k-2.2)**2/0.5)
u_3=1*np.sin(2*np.pi*1*t_k)
# Prepare the signal
Signal=u_1+u_2+u_3+np.random.normal(0,0.1,size=n_t)
plt.plot(t_k,Signal)

# Perform the Frequency Analysis
Signal_FFT = np.fft.fft(Signal)/np.sqrt(n_t) # Compute the DFT
Freqs_S=np.fft.fftfreq(len(Signal))*f_s # Compute the frequency bins



from MRA_Functions_Miguel import MRA_SISO
F_V=np.array([10,70,110,300])
U_MRA,H_MRA=MRA_SISO(Signal, F_V, f_s, 511)


ax1 = plt.subplot(311)
plt.plot(t_k, U_MRA[:,2])
plt.setp(ax1.get_xticklabels(), visible=False)

# share x only
ax2 = plt.subplot(312, sharex=ax1)
plt.plot(t_k, U_MRA[:,3])
# make these tick labels invisible
plt.setp(ax2.get_xticklabels(), visible=False)

# share x and y
ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
plt.plot(t_k, U_MRA[:,4])
plt.xlim(-0.01, 4.1)
plt.show()
plt.savefig('MRA_Results_Signal.pdf', dpi=100)      



# We reconstruct a frequency axis
n_f=np.shape(H_MRA[:,2])[0]
Freqs=np.linspace(0,f_s/2,n_f);



ax1 = plt.subplot(311)
plt.plot(Freqs, H_MRA[:,2])
plt.plot(Freqs_S,abs(Signal_FFT)/max(abs(Signal_FFT)),'k')
plt.setp(ax1.get_xticklabels(), visible=False)

# share x only
ax2 = plt.subplot(312, sharex=ax1)
plt.plot(Freqs, H_MRA[:,3])
plt.plot(Freqs_S,abs(Signal_FFT)/max(abs(Signal_FFT)),'k')

# make these tick labels invisible
plt.setp(ax2.get_xticklabels(), visible=False)

# share x and y
ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
plt.plot(Freqs, H_MRA[:,4])
plt.plot(Freqs_S,abs(Signal_FFT)/max(abs(Signal_FFT)),'k')

plt.xlim(0, 500)
plt.show()
plt.savefig('MRA_Results_Trans_F.pdf', dpi=100)      






