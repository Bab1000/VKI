# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 17:16:44 2020

@author: mendez
"""

# Here I put an home-made function to perform MRA of a 1D signal
import numpy as np
from scipy import signal

def MRA_SISO(u,F_V,f_s,N):
  """
  This function computes the MRA of a signal u
  using Hamming Windows   
  :param u: Input Signal
  :param F_V: Frequency Splitting Vectors (see notes) 
  :param f_s: Sampling Frequency
  :param N: Order of the Filter   
  :return: U_MRA, n_M scale partitions of the signal
  H_MRA  n_M scale Amplitude Responses Functions
  """
  # Get number of scales and number of points
  n_M=len(F_V)+1; n_t=len(u)
  # Initialize the output
  U_MRA=np.zeros((n_t,n_M))
  # Initialize the Transfer Function Matrix
  H_MRA=np.zeros((512,n_M))  
  # Loop from highest frequency to lowest:
  F_V_o=-np.sort(-F_V)
  # Loop over the scales
  for m in range(0,n_M):
   if m==0:
    #Create first low pass kernel     
    h=signal.firwin(N, F_V_o[m], pass_zero=True, fs=f_s)
    # Get frequency response (for checking)
    w,H_L=signal.freqz(h) 
    # This is the first large scale
    u_L=signal.fftconvolve(u,h,'same'); u_H=u-u_L; 
    U_MRA[:,m]=u_H
    H_MRA[:,m]=1-abs(H_L)
    print('Scale '+str(m),'Computed')
   elif m>0 and m<n_M-1:
    #Create mth low pass kernel     
    h=signal.firwin(N, F_V_o[m], pass_zero=True, fs=f_s)
    w,H_L_new=signal.freqz(h) 
    #Low-pass filter m+1
    u_L_new=signal.fftconvolve(u,h,'same')
    # Band Pass filtered and store
    u_H=u_L-u_L_new; U_MRA[:,m]=u_H
    u_L=u_L_new
    # Get Frequency Transf Function and store
    H_MRA[:,m]=abs(H_L)-abs(H_L_new); H_L=H_L_new
    print('Scale '+str(m),'Computed')
   else: 
    U_MRA[:,m]=u_L_new
    # Get Frequency Transf Function and store
    H_MRA[:,m]=abs(H_L); 
    print('Scale '+str(m),'Computed')
    print('Done!')
  return U_MRA,H_MRA