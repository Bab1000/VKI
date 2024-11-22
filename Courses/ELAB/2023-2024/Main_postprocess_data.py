import numpy as np
import os
import csv
#import fig_management
import matplotlib.pyplot as plt
from IR_calibration import IR_calibration
import pandas as pd
from nat_conv_interp import nat_conv_interp

# Specify the path to the folder containing your CSV file
folder_natural_conv = '/home/jpe/VKI/ELAB/2023-2024/Data/Natural_convection'

folder_forced_conv = '/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection'

# Specify the filename of your CSV file
csv_filename1 = 'Natural_conv_17V_1_5A.csv'

# Construct the full path to the CSV file
csv_filepath1 = os.path.join(folder_natural_conv, csv_filename1)

# Load the data from the CSV file using the csv module
with open(csv_filepath1, 'r') as file:
    # Assuming your CSV file has no header and contains numerical values separated by commas
    reader = csv.reader(file, delimiter=';')
    matrix = [[float(num) for num in row] for row in reader]

# Convert the matrix to a NumPy array
image1 = np.array(matrix)

# Specify the filename of your CSV file
csv_filename2 = 'Natural_conv_22V_2_4A.csv'

# Construct the full path to the CSV file
csv_filepath2 = os.path.join(folder_natural_conv, csv_filename2)

# Load the data from the CSV file using the csv module
with open(csv_filepath2, 'r') as file:
    # Assuming your CSV file has no header and contains numerical values separated by commas
    reader2 = csv.reader(file, delimiter=';')
    matrix2 = [[float(num) for num in row] for row in reader2]

# Convert the matrix to a NumPy array
image2 = np.array(matrix2)
# Specify the filename of your CSV file
csv_filename3 = 'Natural_conv_26V_3A.csv'

# Construct the full path to the CSV file
csv_filepath3 = os.path.join(folder_natural_conv, csv_filename3)

# Load the data from the CSV file using the csv module
with open(csv_filepath3, 'r') as file:
    # Assuming your CSV file has no header and contains numerical values separated by commas
    reader3 = csv.reader(file, delimiter=';')
    matrix3 = [[float(num) for num in row] for row in reader3]

# Convert the matrix to a NumPy array
image3 = np.array(matrix3)

# Specify the filename of your CSV file
csv_filename4 = 'Natural_conv_29V_3_2A.csv'

# Construct the full path to the CSV file
csv_filepath4 = os.path.join(folder_natural_conv, csv_filename4)

# Load the data from the CSV file using the csv module
with open(csv_filepath4, 'r') as file:
    # Assuming your CSV file has no header and contains numerical values separated by commas
    reader4 = csv.reader(file, delimiter=';')
    matrix4 = [[float(num) for num in row] for row in reader4]

# Convert the matrix to a NumPy array
image4 = np.array(matrix4)
# Specify the filename of your CSV file
csv_filename5 = 'Natural_conv_34V_4A.csv'

# Construct the full path to the CSV file
csv_filepath5 = os.path.join(folder_natural_conv, csv_filename5)

# Load the data from the CSV file using the csv module
with open(csv_filepath5, 'r') as file:
    # Assuming your CSV file has no header and contains numerical values separated by commas
    reader5 = csv.reader(file, delimiter=';')
    matrix5 = [[float(num) for num in row] for row in reader5]

# Convert the matrix to a NumPy array
image5 = np.array(matrix5)

# data for q balance
V = np.array([17.92, 22.64, 26.4, 29.7, 34.83])
I = np.array([1.5, 2.4, 2.8, 3.2, 3.8])

d = 0.42
eps1 = 0.02
eps2 = 0.95
k = 0.15
t = 1.5e-3
T_amb_nat = 22.6 + 273.15
T_amb_forced = 22.4 + 273.15 
sigma = 5.67e-8
q_joule = 4*V*I/(np.pi*d**2)
q_joule_forced = 120*4/(np.pi*d**2)
rho = 1.225
mu = 1.789e-5
vmin1, vmax1 = 20,70
d_nozzle = 0.0125

"""
im2 = plt.imshow(image, cmap='jet', vmin=vmin1, vmax=vmax1)
plt.title('Original image')
plt.colorbar(im2, label='Values')
"""
n_column = 175
x = np.linspace(0, 320, 321)
y = np.linspace(0, 213, 214)

T_values1 = IR_calibration(image1)[0]
T_values2 = IR_calibration(image2)[0]
T_values3= IR_calibration(image3)[0]
T_values4= IR_calibration(image4)[0]
T_values5= IR_calibration(image5)[0]
"""
im2 = plt.imshow(T_values, cmap='jet', vmin=vmin1, vmax=vmax1)
plt.title('Original image')
plt.colorbar(im2, label='Values')
plt.show()
"""
T_vertical1 = np.flip(T_values1[20:234,n_column] + 273.15,axis=0)
T_vertical2 = np.flip(T_values2[20:234,n_column] + 273.15,axis=0)
T_vertical3 = np.flip(T_values3[20:234,n_column] + 273.15,axis=0)
T_vertical4 = np.flip(T_values4[20:234,n_column] + 273.15,axis=0)
T_vertical5 = np.flip(T_values5[20:234,n_column] + 273.15,axis=0)

print(len(T_vertical1))
print(len(y))

scaling = 0.196


matrix_T = np.column_stack((T_vertical1, T_vertical2, T_vertical3, T_vertical4, T_vertical5))
np.savetxt('matrix_T.csv', matrix_T, delimiter=',', header='T1,T2,T3,T4,T5', comments='')


h1 = (q_joule[0] - (eps1+eps2)*sigma*((T_vertical1)**4-(T_amb_nat)**4))/(2*((T_vertical1)-T_amb_nat))

print(np.shape(q_joule[0] - (eps1+eps2)*sigma*((T_vertical1)**4-(T_amb_nat)**4)))
print(np.shape(2*((T_vertical1)-T_amb_nat)))

h2 = (q_joule[1] - (eps1+eps2)*sigma*((T_vertical2)**4-(T_amb_nat)**4))/(2*((T_vertical2)-T_amb_nat))
h3 = (q_joule[2] - (eps1+eps2)*sigma*((T_vertical3)**4-(T_amb_nat)**4))/(2*((T_vertical3)-T_amb_nat))
h4 = (q_joule[3] - (eps1+eps2)*sigma*((T_vertical4)**4-(T_amb_nat)**4))/(2*((T_vertical4)-T_amb_nat))
h5 = (q_joule[4] - (eps1+eps2)*sigma*((T_vertical5)**4-(T_amb_nat)**4))/(2*((T_vertical5)-T_amb_nat))

matrix_h = np.column_stack((h1,h2,h3,h4,h5))
np.savetxt('matrix_h.csv', matrix_h, delimiter=',', header='h1,h2,h3,h4,h5', comments='')

plt.figure()
plt.plot(T_vertical1,y*scaling, label = '$q_{\\mathrm{J},1}$' )
plt.plot(T_vertical2,y*scaling, label='$q_{\\mathrm{J},2}$' )
plt.plot(T_vertical3,y*scaling, label='$q_{\\mathrm{J},3}$')
plt.plot(T_vertical4,y*scaling, label = '$q_{\\mathrm{J},4}$')
plt.plot(T_vertical5,y*scaling, label = '$q_{\\mathrm{J},5}$', linewidth=1)
plt.xlabel("$T$ [K]")
plt.ylabel("$y$ [cm]")
plt.xlim(left=298)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig("/home/jpe/VKI/ELAB/2023-2024/Plots/natural_T.pdf")

plt.figure()
plt.plot(h1,y*scaling , label='$q_\\mathrm{J} = ' +str(round(q_joule[0], 2)) + '~\\mathrm{W/m^2}$')
plt.plot(h2,y*scaling, label='$q_\\mathrm{J} = ' +str(round(q_joule[1], 2)) + '~\\mathrm{W/m^2}$')
plt.plot(h3,y*scaling, label='$q_\\mathrm{J} = ' +str(round(q_joule[2], 2)) + '~\\mathrm{W/m^2}$')
plt.plot(h4,y*scaling, label='$q_\\mathrm{J} = ' +str(round(q_joule[3], 2)) + '~\\mathrm{W/m^2}$')
plt.plot(h5,y*scaling, label='$q_\\mathrm{J} = ' +str(round(q_joule[4], 2)) + '~\\mathrm{W/m^2}$', linewidth=1)
plt.xlabel("$h$ [W/m$^2\\cdot$K]")
plt.ylabel("$y$ [cm]")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.savefig("/home/jpe/VKI/ELAB/2023-2024/Plots/natural_h.pdf")

# np.flip(matrix, axis=0) to fli matrix 
forced1 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_1D_6.csv", sep=';').to_numpy()
forced2 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_1D_12.csv", sep=';').to_numpy()
forced3 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_1D_18.csv", sep=';').to_numpy()
forced4 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_5D_6.csv", sep=';').to_numpy()
forced5 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_5D_12.csv", sep=';').to_numpy()
forced6 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_5D_18.csv", sep=';').to_numpy()
forced7 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_5D_18.csv", sep=';').to_numpy()
forced8 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_5D_18.csv", sep=';').to_numpy()
forced9 = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/Data/Forced_convection/Forced_conv_5D_18.csv", sep=';').to_numpy()

forced1 = forced1.astype(float)
forced2 = forced2.astype(float)
forced3 = forced3.astype(float)
forced4 = forced4.astype(float)
forced5 = forced5.astype(float)
forced6 = forced6.astype(float)
forced7 = forced7.astype(float)
forced8 = forced8.astype(float)
forced9 = forced9.astype(float)

T_forced1 = IR_calibration(forced1)[0]
print(T_forced1[0][0])
T_forced2 = IR_calibration(forced2)[0]
T_forced3 = IR_calibration(forced3)[0]
T_forced4 = IR_calibration(forced4)[0]
T_forced5 = IR_calibration(forced5)[0]
T_forced6 = IR_calibration(forced6)[0]
T_forced7 = IR_calibration(forced7)[0]
T_forced8 = IR_calibration(forced8)[0]
T_forced9 = IR_calibration(forced9)[0]


T_vertical_forced1 = np.flip(T_forced1[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced2 = np.flip(T_forced2[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced3 = np.flip(T_forced3[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced4 = np.flip(T_forced4[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced5 = np.flip(T_forced5[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced6 = np.flip(T_forced6[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced7 = np.flip(T_forced7[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced8 = np.flip(T_forced8[20:234,n_column] + 273.15, axis = 0)
T_vertical_forced9 = np.flip(T_forced9[20:234,n_column] + 273.15, axis = 0)



#print(min(T_vertical_forced3), max(T_vertical_forced3))

# compute h_natural for each y and measurement

# Placeholder for h_nat values
h_nat_1= np.zeros_like(y)
h_nat_2= np.zeros_like(y)
h_nat_3= np.zeros_like(y)
h_nat_4= np.zeros_like(y)
h_nat_5= np.zeros_like(y)
h_nat_6= np.zeros_like(y)
h_nat_7= np.zeros_like(y)
h_nat_8= np.zeros_like(y)
h_nat_9= np.zeros_like(y)

for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced1[0:i+1]))  # Assuming T is a float value
    h_nat_1[i] = nat_conv_interp(y_index, T_value)





for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced2[0:i+1]))  # Assuming T is a float value
    h_nat_2[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced3[0:i+1]))  # Assuming T is a float value
    h_nat_3[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced4[0:i+1]))  # Assuming T is a float value
    h_nat_4[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced5[0:i+1]))  # Assuming T is a float value
    h_nat_5[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced6[0:i+1]))  # Assuming T is a float value
    h_nat_6[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced7[0:i+1]))  # Assuming T is a float value
    h_nat_7[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced8[0:i+1]))  # Assuming T is a float value
    h_nat_8[i] = nat_conv_interp(y_index, T_value)



for i in range(len(y)):
    y_index = int(y[i])  # Assuming y is an index (integer)
    T_value = float(np.mean(T_vertical_forced9[0:i+1]))  # Assuming T is a float value
    h_nat_9[i] = nat_conv_interp(y_index, T_value)


T_for = [T_vertical_forced1, T_vertical_forced2, T_vertical_forced3, T_vertical_forced4, T_vertical_forced5, T_vertical_forced6, T_vertical_forced7, T_vertical_forced8, T_vertical_forced9]
h_for = [h_nat_1, h_nat_2, h_nat_3,h_nat_4, h_nat_5, h_nat_6,h_nat_7, h_nat_8, h_nat_9]

T_spiral = []
h_forced = []

# heat balance computations forced convection
for j in range(len(T_for)):
    T_spiral.append((t/k)*(eps2*sigma*(T_for[j]**4 - T_amb_forced**4)+h_for[j]*(T_for[j]-T_amb_forced)) + T_for[j])
    h_forced.append((q_joule_forced - eps1*sigma*(T_spiral[j]**4-T_amb_forced**4) - eps2*sigma*(T_for[j]**4-T_amb_forced**4)-h_for[j]*(T_for[j]-T_amb_forced))/(T_spiral[j]-T_amb_forced))

"""
plt.figure()
plt.plot(h_nat_6,y*scaling)
plt.xlabel("natural forced $h$ [W/m$^2$K]")
plt.ylabel("$y$ [cm]")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig("JAN_natural_forced_h.pdf")

plt.figure()
plt.plot(T_spiral[5],y*scaling)
plt.xlabel("Spiral temperature [K]")
plt.ylabel("$y$ [cm]")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig("JAN_natural_forced_T1.pdf")

plt.figure()
plt.plot(T_for[5],y*scaling)
plt.xlabel("Measured temperature [K]")
plt.ylabel("$y$ [cm]")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig("JAN_natural_forced_T2.pdf")
"""

np.savetxt("matT_spiral.csv", np.array(T_spiral), delimiter=',')
np.savetxt("matT_for.csv", np.array(T_for), delimiter=',')

x_rotameter = [60,120,180]
Re = []
v_flow_rate = []
for jj in range(len(x_rotameter)):
    v_flow_rate.append((0.0129*x_rotameter[jj] + 0.0308)/1000)
    Re.append(4*v_flow_rate[jj]*rho/(np.pi*d_nozzle*mu))


plt.figure()
plt.plot(h_forced[0],y*scaling , label='$Re = ' + str(round(Re[0],2)) + '$')
plt.plot(h_forced[1],y*scaling , label='$Re = ' + str(round(Re[1],2)) + '$')
plt.plot(h_forced[2],y*scaling , label='$Re = ' + str(round(Re[2],2)) + '$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$h$ [W/(m$^2\\cdot$K)]")
plt.ylabel("$y$ [cm]")
plt.savefig("/home/jpe/VKI/ELAB/2023-2024//Plots/forced_h_dist1.pdf")

plt.figure()
plt.plot(h_forced[3],y*scaling , label='$Re = ' + str(round(Re[0],2)) + '$')
plt.plot(h_forced[4],y*scaling , label='$Re = ' + str(round(Re[1],2)) + '$')
plt.plot(h_forced[5],y*scaling , label='$Re = ' + str(round(Re[2],2)) + '$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$h$ [W/(m$^2\\cdot$K)]")
plt.ylabel("$y$ [cm]")
plt.savefig("/home/jpe/VKI/ELAB/2023-2024//Plots/forced_h_dist2.pdf")

plt.figure()
plt.plot(h_forced[6],y*scaling , label='$Re = ' + str(round(Re[0],2)) + '$')
plt.plot(h_forced[7],y*scaling , label='$Re = ' + str(round(Re[1],2)) + '$')
plt.plot(h_forced[8],y*scaling , label='$Re = ' + str(round(Re[2],2)) + '$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$h$ [W/(m$^2\\cdot$K)]")
plt.ylabel("$y$ [cm]")
plt.savefig("/home/jpe/VKI/ELAB/2023-2024//Plots/forced_h_dist3.pdf")

plt.figure()
plt.plot(T_vertical_forced1,y*scaling , label='$Re = ' + str(round(Re[0],2)) + '$')
plt.plot(T_vertical_forced2,y*scaling , label='$Re = ' + str(round(Re[1],2)) + '$')
plt.plot(T_vertical_forced3,y*scaling , label='$Re = ' + str(round(Re[2],2)) + '$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$T$ [K]")
plt.ylabel("$y$ [cm]")
plt.savefig("/home/jpe/VKI/ELAB/2023-2024//Plots/forced_T_dist1.pdf")

plt.figure()
plt.plot(T_vertical_forced4,y*scaling , label='$Re = ' + str(round(Re[0],2)) + '$')
plt.plot(T_vertical_forced5,y*scaling , label='$Re = ' + str(round(Re[1],2)) + '$')
plt.plot(T_vertical_forced6,y*scaling , label='$Re = ' + str(round(Re[2],2)) + '$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$T$ [K]")
plt.ylabel("$y$ [cm]")
plt.savefig("/home/jpe/VKI/ELAB/2023-2024//Plots/forced_T_dist2.pdf")

plt.figure()
plt.plot(T_vertical_forced7,y*scaling , label='$Re = ' + str(round(Re[0],2)) + '$')
plt.plot(T_vertical_forced8,y*scaling , label='$Re = ' + str(round(Re[1],2)) + '$')
plt.plot(T_vertical_forced9,y*scaling , label='$Re = ' + str(round(Re[2],2)) + '$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("$T$ [K]")
plt.ylabel("$y$ [cm]")
plt.savefig("/home/jpe/VKI/ELAB/2023-2024//Plots/forced_T_dist3.pdf")