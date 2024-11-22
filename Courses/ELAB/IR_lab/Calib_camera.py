import numpy as np
from scipy.optimize import least_squares
#import fig_management
import matplotlib.pyplot as plt

#plt.style.use("grayscale")

# Define the correlation equation
def correlation_equation(coefficients, T):
    R, B, F = coefficients
    return R / (np.exp(B / T) - F)

# Thermocouples data 
TC1_data = np.array([71.64333333, 65.16, 56.77666667, 42.34333333, 36.54, 29.96666667])
TC2_data = np.array([69.82666667, 62.63333333, 55.7, 41.68, 35.99666667, 29.53666667])

# TC calibration coefficients (from the linear fit)
# new
# a1 = 1.0094 
# a2 = 1.0090
# b1 = -4.4428
# b2 = -4.4854

# old
a1 = 0.9999 
a2 = 0.9981
b1 = -2.6568
b2 = -2.5917

# Computing the real temperature from thermocouples
T1_real_data =  a1*TC1_data + b1 
T2_real_data = a2*TC2_data + b2

# IU data from the camera calibration 
IU_data_TC1 = np.array([ 8516.333333, 7704, 6781.666667, 5282, 4750, 4184.333333])  
IU_data_TC2 = np.array([ 8598.666667, 7707.666667, 6766.333333, 5270.333333, 4767.666667, 4170.666667]) 


# Initial coefficients for the least square fit
initial_coefficients = [-2000, -40, 1]

# Define the objective function for least squares
def objective_function(coefficients):
    return correlation_equation(coefficients, T1_real_data) - IU_data_TC1

# Fit the correlation equation to the data using least squares
result1 = least_squares(objective_function, initial_coefficients)

# Define the objective function for least squares
def objective_function(coefficients):
    return correlation_equation(coefficients, T2_real_data) - IU_data_TC2

# Fit the correlation equation to the data using least squares
result2 = least_squares(objective_function, initial_coefficients)


# Extract the optimized coefficients
R1_optimized, B1_optimized, F1_optimized = result1.x     # .x is the coefficients of the best fitting
R2_optimized, B2_optimized, F2_optimized = result2.x
# Print the optimized coefficients
print("Optimized R1, R2:", R1_optimized, R2_optimized )
print("Optimized B1, B2:", B1_optimized, B2_optimized)
print("Optimized F1, F2:", F1_optimized, F2_optimized)

IU_trial1 = R1_optimized/(np.exp(B1_optimized/55.2) - F1_optimized)
IU_trial2 = R2_optimized/(np.exp(B2_optimized/54.5) - F2_optimized)

print("IU_trial:", IU_trial1, IU_trial2)

# just for the plot (x-axis)
x_T = np.arange(20,75)


plt.figure()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.plot(x_T, ((R1_optimized / (np.exp(B1_optimized / x_T) - F1_optimized)) +
               (R2_optimized / (np.exp(B2_optimized / x_T) - F2_optimized))) / 2,
         linestyle=":", color="black", label='Average calibration')
plt.plot(x_T, R1_optimized / (np.exp(B1_optimized / x_T) - F1_optimized),
         linestyle="-.", color="black", label='Calibration from TC 1')
plt.plot(x_T, R2_optimized / (np.exp(B2_optimized / x_T) - F2_optimized),
         linestyle="-", color="black", label='Calibration from TC 2')
plt.scatter(T1_real_data, IU_data_TC1, marker='o', color="gray", s=20, label='TC 1 Temperatures')  # 's' controls the marker size
plt.scatter(T2_real_data, IU_data_TC2, marker='^', color="gray", s=20, label='TC 2 temperatures')
plt.xlabel("$T$ [$^{\\circ}$C]")
plt.ylabel("Camera counts [IU]")
plt.legend()
plt.savefig("/home/jpe/VKI/ELAB/IR_lab/Plots/Camera_calib.pdf")
plt.close()




