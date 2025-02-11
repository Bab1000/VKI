import numpy as np
from Utils_SP_Exercises_Mendez import *
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)

# Path to the data
path = "/home/jpe/VKI/Courses/SP/Mendez_part"

# --------------------------------------------------------------------------

# Check for the libraries
# -----------------------

LibraryManager()

# --------------------------------------------------------------------------

# Exercise 1
# ----------
print("")
print(Fore.BLUE + "Question 1 :")
print(Fore.BLUE + "------------")

# Inputs
fs = 3000   # Sampling frequency
n_samples = 13200   # Number of samples
coord = [19.444,9.0254]   # Coordinate at which the spectral analysis si required

# Processing the data
u,v = ResFileProcessing(path,coord)

# Exercise 1
freq, fft, vel_magnitude,vel_magnitude_filtered = Exercise1(u,v,fs,n_samples)

# Plots
PlotExercise1(freq,fft,path)
print("")

# --------------------------------------------------------------------------

# Exercise 2
# ----------

print("")
print(Fore.BLUE + "Question 2 :")
print(Fore.BLUE + "------------")

freq, t, Amplitude = Exercise2(vel_magnitude_filtered,fs)

# Plots
PlotExercise2(freq,t,Amplitude,path)

print("")

# --------------------------------------------------------------------------

# Exercise 3
# ----------

print("")
print(Fore.BLUE + "Question 3 :")
print(Fore.BLUE + "------------")

U_MRA, H_MRA, freq_splitting_vector= Exercise3(vel_magnitude,fs,n_samples)

# Plots
PlotExercise3(U_MRA,freq_splitting_vector,n_samples,fs,path)
print("")
print(Fore.BLUE + "End of the script ! ")
print("")