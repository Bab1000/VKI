import numpy as np
#import fig_management
import matplotlib.pyplot as plt
from IR_calibration import IR_calibration
import pandas as pd
from nat_conv_interp import nat_conv_interp

# INPUTS 
# ------

# diameter of the heated surface
d = 0.46

# emissivity of the spiral
eps1 = 0.02

# epoxy plate emissivity
eps2 = 0.95

# thermal conductivity
k = 0.15

# epoxy plate thickness
t = 1.5e-3

# ambiant temperature for NATURAL convection experimental campaign
T_amb_nat = 22.6 + 273.15

# ambiant temperature for FORCED convection experimental campaign
T_amb_forced = 22.4 + 273.15

# Stephan-Boltzmann cst
sigma = 5.67e-8

#heat flux due to Joul effect
P = 120     # max power without exceeding the max current
n = 4       # number of spirals ??
q_joule_forced = P*n/(np.pi*d**2)

# density of air
rho = 1.225

#dynamic viscosity of air
mu = 1.789e-5

# nbr of collumns in the csv files
n_column = 175

# number of columns in csv files for forced convection    ( linked to the window utilized later for calibration ln: 95-101 )
x = np.linspace(0, 320, 321)
# number of lines in csv files for forced convection      ( linked to the window utilized later for calibration ln: 95-101 )
y = np.linspace(0, 213, 214)

# scaling for plot
scaling = 0.196

# CSV DATA DIRECTORY + READING
# ----------------------------

# reading data
forced_5 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_5deg.csv", sep=';').to_numpy()
#forced_10 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_10deg.csv", sep=';').to_numpy()
forced_15 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_15deg.csv", sep=';').to_numpy()
forced_neg5 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_-5deg.csv", sep=';').to_numpy()
#forced_neg10 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_-10deg.csv", sep=';').to_numpy()
forced_neg15 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_-15deg.csv", sep=';').to_numpy()
forced_0 = pd.read_csv("/home/jpe/VKI/ELAB/Data/Forced_convection/Forced_conv_5D_12.csv", sep=';').to_numpy()

# ===========================================================================================================
# ===========================================================================================================
# ==========================================================================================================

# Converting integers into floats
# -------------------------------

forced_5 = forced_5.astype(float)
#forced_10 = forced_10.astype(float)
forced_15 = forced_15.astype(float)
forced_neg5 = forced_neg5.astype(float)
#forced_neg10 = forced_neg10.astype(float)
forced_neg15 = forced_neg15.astype(float)
forced_0 = forced_0.astype(float)

# Computing temperature from data
# --------------------------------

T_forced_5 = IR_calibration(forced_5)[0]  # [0] means that only the first element of the vectors is returned which is T_avg
#T_forced_10 = IR_calibration(forced_10)[0]
T_forced_15 = IR_calibration(forced_15)[0]
T_forced_neg5 = IR_calibration(forced_neg5)[0]
#T_forced_neg10 = IR_calibration(forced_neg10)[0]
T_forced_neg15 = IR_calibration(forced_neg15)[0]
T_forced_0 = IR_calibration(forced_0)[0]

T_vertical_forced_5 = np.flip(T_forced_5[20:234,n_column] + 273.15, axis = 0)               # fliping because the data is from top of image to bottom and we want opposite
T_vertical_forced_10 = 0 #np.flip(T_forced_10[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced_15 = np.flip(T_forced_15[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced_neg5 = np.flip(T_forced_neg5[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced_neg10 = 0 #np.flip(T_forced_neg10[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced_neg15 = np.flip(T_forced_neg15[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced_0 = np.flip(T_forced_0[20:234,n_column] + 273.15, axis = 0)

# compute h_natural for each y and measurement

# Placeholder for h_nat values
h_nat_1= np.zeros_like(y)
h_nat_2= np.zeros_like(y)
h_nat_3= np.zeros_like(y)
h_nat_4= np.zeros_like(y)
h_nat_5= np.zeros_like(y)
h_nat_6= np.zeros_like(y)
h_nat_7= np.zeros_like(y)


for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced_5[0:i+1]))  # Assuming T is a float value
    h_nat_1[i] = nat_conv_interp(y_index, T_value)
   

#for i in range(len(y)):
    #y_index = int(y[i])  # Assuming y is an index (integer)
    #T_value = float(np.mean(T_vertical_forced_10[0:i+1]))  # Assuming T is a float value
    #h_nat_2[i] = nat_conv_interp(y_index, T_value)

for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced_15[0:i+1]))  # Assuming T is a float value
    h_nat_3[i] = nat_conv_interp(y_index, T_value)

for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced_neg5[0:i+1]))  # Assuming T is a float value
    h_nat_4[i] = nat_conv_interp(y_index, T_value)

#for i in range(len(y)):
    #y_index = int(y[i])  # Assuming y is an index (integer)
    #T_value = float(np.mean(T_vertical_forced_neg10[0:i+1]))  # Assuming T is a float value
    #h_nat_5[i] = nat_conv_interp(y_index, T_value)

for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced_neg15[0:i+1]))  # Assuming T is a float value
    h_nat_6[i] = nat_conv_interp(y_index, T_value)

for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced_0[0:i+1]))  # Assuming T is a float value
    h_nat_7[i] = nat_conv_interp(y_index, T_value)



T_for = [T_vertical_forced_5, T_vertical_forced_10, T_vertical_forced_15, T_vertical_forced_neg5, T_vertical_forced_neg10, T_vertical_forced_neg15, T_vertical_forced_0]
h_for = [h_nat_1, h_nat_2, h_nat_3,h_nat_4, h_nat_5, h_nat_6,h_nat_7]

T_spiral = []
h_forced = []

# heat balance computations forced convection
for j in range(len(T_for)):
    T_spiral.append((t/k)*(eps2*sigma*(T_for[j]**4 - T_amb_forced**4)+h_for[j]*(T_for[j]-T_amb_forced)) + T_for[j])
    h_forced.append((q_joule_forced - eps1*sigma*(T_spiral[j]**4-T_amb_forced**4) - eps2*sigma*(T_for[j]**4-T_amb_forced**4)-h_for[j]*(T_for[j]-T_amb_forced))/(T_spiral[j]-T_amb_forced))

plt.figure()
plt.plot(h_forced[2],y*scaling,label='$\\alpha = 15^{\\circ}$')
plt.plot(h_forced[0],y*scaling, label='$\\alpha = 5^{\\circ}$')
plt.plot(h_forced[6],y*scaling,label='$\\alpha = 0^{\\circ}$')
plt.plot(h_forced[5],y*scaling,label='$\\alpha = -15^{\\circ}$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$h$ [W/(m$^2\\cdot$K)]")
plt.ylabel('$y$ [cm]')
plt.savefig("/home/jpe/VKI/ELAB/Plots/Plotsforced_angles_h.pdf")