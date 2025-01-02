import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from CSV_reader import CSV_reader
import pdb
from h_correlation import *

plt.style.use("grayscale")

# ========================================================================

# READING NATURAL/FORCED CONVECTION FILES
# ---------------------------------------

# Specify the path to the folder containing your CSV file
folder_natural_conv = '/home/jpe/VKI/Courses/ELAB/IR_lab/Data/Natural_convection'

Images_natural = {}

Images_natural = CSV_reader(Images_natural,folder_natural_conv)

# ========================================================================

# INPUTS
# ------

# data for q balance
V = np.array([10.5, 15.2, 19.9, 25.1, 28.5])
I = np.array([1.16, 1.67, 2.14, 2.64, 2.94])

d = 0.46
eps1 = 0.02
eps2 = 0.95
k = 0.15
t = 1.5e-3
T_amb_nat = 19.6 + 273.15
sigma = 5.67e-8
q_joule = 4*V*I/(np.pi*d**2)
rho = 1.225
mu = 1.789e-5
vmin1, vmax1 = 20,70

# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

# ---------------
# | COMPUTATION |
# ---------------

# ========================================================================

# Natural convection
# ------------------

T_list_naturalV = []
keys_list_naturalV = []
h_list_naturalV = []
T_list_naturalH = []
keys_list_naturalH = []
h_list_naturalH = []
count_h_natural = 0

# with new calibration
# R1, R2 = -2822.9389090622676, -2754.82002517215
# B1, B2 = -39.58633426427951, -41.676182399554875
# F1, F2 = 0.8899176908949745, 0.85144000173110

# with old calib
R1, R2 = -2705.1641333102857, -2640.6988306265257
B1, B2 = -38.97839012198067, -41.20604904144822
F1, F2 = 0.8862613818859224, 0.8476297258370964

for key in sorted(Images_natural.keys()):

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
    
    if "vertical" in key:
        
        # Listing vertical temperatures for the differnt images
        T_list_naturalV.append(T_avgT)
        keys_list_naturalV.append(key)

        # compute h
        h_naturalV = (q_joule[count_h_natural] - (eps1+eps2)*sigma*((T_avgh)**4-(T_amb_nat)**4))/(2*((T_avgh)-T_amb_nat))

        # listing h
        h_list_naturalV.append(h_naturalV)
        count_h_natural = count_h_natural + 1
    else:
        # Listing vertical temperatures for the differnt images
        T_list_naturalH.append(T_avgT)
        keys_list_naturalH.append(key)

        # compute h
        h_naturalH = (q_joule[count_h_natural] - (eps1+eps2)*sigma*((T_avgh)**4-(T_amb_nat)**4))/(2*((T_avgh)-T_amb_nat))
    
        # listing h
        h_list_naturalH.append(h_naturalH)
        
        count_h_natural = count_h_natural + 1

    if (count_h_natural == 5):
        count_h_natural = 0

# create the T matrix csv file
T_list_naturalV = np.array(T_list_naturalV)
T_list_naturalH = np.array(T_list_naturalH)
matrix_TV = np.column_stack(T_list_naturalV)
matrix_TH = np.column_stack(T_list_naturalH)
np.savetxt('/home/jpe/VKI/Courses/ELAB/IR_lab/Matrices/matrix_TV.csv', matrix_TV, delimiter=',', header=','.join(keys_list_naturalV), comments='')
np.savetxt('/home/jpe/VKI/Courses/ELAB/IR_lab/Matrices/matrix_TH.csv', matrix_TH, delimiter=',', header=','.join(keys_list_naturalH), comments='')

# create the h matrix csv file
h_list_naturalV = np.array(h_list_naturalV)
h_list_naturalH = np.array(h_list_naturalH)
matrix_hV = np.column_stack(h_list_naturalV)
matrix_hH = np.column_stack(h_list_naturalH)
np.savetxt('/home/jpe/VKI/Courses/ELAB/IR_lab/Matrices/matrix_hV.csv', matrix_hV, delimiter=',', header=','.join(keys_list_naturalV), comments='')
np.savetxt('/home/jpe/VKI/Courses/ELAB/IR_lab/Matrices/matrix_hH.csv', matrix_hH, delimiter=',', header=','.join(keys_list_naturalH), comments='')

# ========================================================================

# Literature comparison
# ---------------------

T_full = np.array(matrix_TV[:,3])
T_cut = T_full[0:120]
T_avg = np.mean(T_cut)

h_natural_lit = h_natural(T_avg)

h_full = np.array(matrix_hV[:,3])
h_cut = h_full[0:120]
h_avg = np.mean(h_cut)

err = np.abs(1 - h_avg/h_natural_lit)*100

print()
print("Natural convecion:")
print("------------------")
print()
print("T : " + str(T_avg))
print("h literature : " + str(h_natural_lit))
print("h computed : " + str(h_avg))
print("Relative error : " + str(err) + " %")

T_forced_mean = 295.46
vdot = 0.00236
dist = 0.110

h_forced_lit = h_forced(T_forced_mean,vdot,dist)

h_forced_comp = 138.69

err = np.abs(1 - h_forced_comp/h_forced_lit)*100

print()
print("Forced convecion:")
print("------------------")
print()
print("T : " + str(T_forced_mean))
print("h literature : " + str(h_forced_lit))
print("h computed : " + str(h_forced_comp))
print("Relative error : " + str(err) + " %")


# ========================================================================
# ========================================================================
# ========================================================================
# ========================================================================

# ---------
# | Plots |
# ---------


# ========================================================================

# scaling for plots
scaling = 2

yTV = np.linspace(0, matrix_TV.shape[0] - 1, matrix_TV.shape[0])
yhV = np.linspace(0, matrix_hV.shape[0] - 1, matrix_hV.shape[0])
yTH = np.linspace(0, matrix_TH.shape[0] - 1, matrix_TH.shape[0])
yhH = np.linspace(0, matrix_hH.shape[0] - 1, matrix_hH.shape[0])

t_amb = np.ones((matrix_TV.shape[0])) * T_amb_nat

plt.figure()
for i in range(matrix_TV.shape[1]):  # Loop over columns
    TV = matrix_TV[:, i]  # Access the i-th column
    if (i != 0):
        plt.plot(TV, yTV * scaling, label='$q_\\mathrm{J} = ' + str(round(q_joule[i], 2)) + '~\\mathrm{W/m^2}$')

#plt.plot(t_amb, yTV * scaling)
plt.xlabel("$T$ [K]")
plt.ylabel("$y$ [mm]")
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig("/home/jpe/VKI/Courses/ELAB/IR_lab/Plots/natural_TV.pdf")

plt.figure()
for i in range(matrix_hV.shape[1]):  # Loop over columns
    hV = matrix_hV[:, i]  # Access the i-th column
    if (i != 0):
        plt.plot(hV, yhV * scaling, label='$q_\\mathrm{J} = ' + str(round(q_joule[i], 2)) + '~\\mathrm{W/m^2}$')

plt.xlabel("$h$ [W/m$^2\\cdot$K]")
plt.ylabel("$y$ [mm]")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.savefig("/home/jpe/VKI/Courses/ELAB/IR_lab/Plots/natural_hV.pdf")

plt.figure()
for i in range(matrix_TH.shape[1]):  # Loop over columns
    TH = matrix_TH[:, i]  # Access the i-th column
    plt.plot(TH, yTH * scaling, label='$q_\\mathrm{J} = ' + str(round(q_joule[i], 2)) + '~\\mathrm{W/m^2}$')


plt.xlabel("$T$ [K]")
plt.ylabel("$x$ [mm]")
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig("/home/jpe/VKI/Courses/ELAB/IR_lab/Plots/natural_TH.pdf")
plt.figure()

plt.figure()
for i in range(matrix_hH.shape[1]):  # Loop over columns
    hH = matrix_hH[:, i]  # Access the i-th column
    plt.plot(hH, yhH * scaling, label='$q_\\mathrm{J} = ' + str(round(q_joule[i], 2)) + '~\\mathrm{W/m^2}$')

plt.xlabel("$h$ [W/m$^2\\cdot$K]")
plt.ylabel("$x$ [mm]")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.savefig("/home/jpe/VKI/Courses/ELAB/IR_lab/Plots/natural_hH.pdf")