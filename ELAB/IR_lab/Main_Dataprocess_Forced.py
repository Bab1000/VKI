import numpy as np
#import fig_management
import matplotlib.pyplot as plt
from IR_calibration import IR_calibration
import pandas as pd
#from nat_conv_interp import nat_conv_interp
from CSV_reader import CSV_reader
from nat_conv_interp import nat_conv_interp
import pdb

# ========================================================================

# READING NATURAL/FORCED CONVECTION FILES
# ---------------------------------------

# Specify the path to the folder containing your CSV file
folder_natural_conv = '/home/jpe/VKI/ELAB/IR_lab/Data/Natural_convection'

folder_forced_conv = '/home/jpe/VKI/ELAB/IR_lab/Data/Forced_convection'

Images_natural = {}
Images_forced = {}

Images_natural = CSV_reader(Images_natural,folder_natural_conv)
Images_forced = CSV_reader(Images_forced,folder_forced_conv)

# ========================================================================

# INPUTS
# ------

# data for q balance
V = np.array([17.92, 22.64, 26.4, 29.7, 34.83])
I = np.array([1.5, 2.4, 2.8, 3.2, 3.8])

d = 0.42
eps1 = 0.02
eps2 = 0.95
k = 0.15
t = 1.5e-3
T_amb_nat = 19.6 + 273.15
T_amb_forced = 19.6 + 273.15 
sigma = 5.67e-8
q_joule = 4*V*I/(np.pi*d**2)
q_joule_forced = 120*4/(np.pi*d**2)
rho = 1.225
mu = 1.789e-5
vmin1, vmax1 = 20,70
d_nozzle = 0.0125

n_column = 175
x = np.linspace(0, 320, 321)
y = np.linspace(0, 213, 214)

# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

# ---------------
# | COMPUTATION |
# ---------------


# ========================================================================

# Forced convection
# ------------------

T_list_forced = []
keys_list_forced = []
h_list_forced = []
count_h_nat = 0

# with old calib
R1, R2 = -2705.1641333102857, -2640.6988306265257
B1, B2 = -38.97839012198067, -41.20604904144822
F1, F2 = 0.8862613818859224, 0.8476297258370964

for key in sorted(Images_forced.keys()):

     # gathering the IU
    line = Images_natural[key]

    # Convert the array to strings first
    line = line.astype(str)
    # Now replace commas with dots
    line = np.char.replace(line, ',', '.')
    # Convert the array to float after the replacement
    line = np.array(line).astype(float)

    # gathering temperature
    T_avgT = ((B1 / np.log((R1/line[0:207]) + F1)) + (B2 / np.log((R2/line[0:207]) + F2)))/2 + 273.15
    T_avgh = ((B1 / np.log((R1/line[:200]) + F1)) + (B2 / np.log((R2/line[:200]) + F2)))/2 + 273.15

    # vertical temperature
    T_vertical_forced = T_avgT + 273.15
    print(T_vertical_forced[0])

    # Listing vertical temperatures for the differnt images
    T_list_forced.append(T_vertical_forced)
    keys_list_forced.append(key)

    h_nat = [0] * len(y)  # array of len y with zeros

    # compute h
    for i in range(len(y)):
        y_index = int(y[i])  # Assuming y is an index (integer)
        T_value = float(np.mean(T_vertical_forced[0:i+1]))  # Assuming T is a float value
        h_nat[i] = nat_conv_interp(y_index, T_value)

    # Listing vertical temperatures for the differnt imagess
    h_list_forced.append(h_nat)

T_list_forced = np.array(T_list_forced)
print(T_list_forced.shape)  # This will now work

matrix_T_forced = np.column_stack(T_list_forced)
print(matrix_T_forced.shape)

T_spiral = []
h_forced = []

# heat balance computations forced convection
for j in range(matrix_T_forced.shape[1]):
    T_spiral.append((t/k)*(eps2*sigma*(matrix_T_forced[:,j]**4 - T_amb_forced**4)+h_list_forced[j]*(matrix_T_forced[:,j]-T_amb_forced)) + matrix_T_forced[:,j])
    h_forced.append((q_joule_forced - eps1*sigma*(T_spiral[j]**4-T_amb_forced**4) - eps2*sigma*(matrix_T_forced[:,j]**4-T_amb_forced**4)-h_list_forced[j]*(matrix_T_forced[:,j]-T_amb_forced))/(T_spiral[j]-T_amb_forced))

# save matrices
np.savetxt("matT_spiral_new.csv", np.array(T_spiral), delimiter=',')
np.savetxt("matT_forced_new.csv", matrix_T_forced, delimiter=',')