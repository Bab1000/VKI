import numpy as np
from scipy.optimize import least_squares
#import fig_management
import matplotlib.pyplot as plt

def IR_calibration(IU):   # IU = intensity units returned by the camera

    # Define the correlation equation
    def correlation_equation(coefficients, T):
        R, B, F = coefficients
        return R / (np.exp(B / T) - F)

    # Thermocouples data ( !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! )
    TC1_data = np.array([68.6, 64.2, 59.8, 55.3, 50.3, 45.3, 40.2, 34.7, 29.7, 23.8])
    TC2_data = np.array([67.3, 63, 59.1, 54.5, 49.8, 44.8, 40, 34.6, 29.8, 23.8])

    # TC calibration coefficients (from the linear fit)
    a1 = 1.008 
    a2 = 1.007
    b1 = 0.004
    b2 = 0.005

    # Computing the real temperature from thermocouples
    T1_real_data =  a1*TC1_data + b1 
    T2_real_data = a2*TC2_data + b2

    # IU data from the camera calibration ( !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ) 
    IU_data = np.array([ 8947, 8338, 7757, 7150, 6555, 5990, 5446, 4895, 4422, 3894])  # Replace with your corresponding A data

    # Initial coefficients for the least square fit
    initial_coefficients = [-2000.0, -40.0, 1.0]

    # Define the objective function for least squares
    def objective_function(coefficients):
        return correlation_equation(coefficients, T1_real_data) - IU_data

    # Fit the correlation equation to the data using least squares
    result1 = least_squares(objective_function, initial_coefficients)

    # Define the objective function for least squares
    def objective_function(coefficients):
        return correlation_equation(coefficients, T2_real_data) - IU_data

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
    plt.plot(x_T,((R1_optimized/(np.exp(B1_optimized/x_T) - F1_optimized))+(R2_optimized/(np.exp(B2_optimized/x_T) - F2_optimized)))/2, label='Avg calibration')
    plt.plot(x_T,R1_optimized/(np.exp(B1_optimized/x_T) - F1_optimized), label='Calibration fit from TC 1')
    plt.plot(x_T,R2_optimized/(np.exp(B2_optimized/x_T) - F2_optimized), label='Calibration fit from TC 2')
    plt.scatter(TC1_data, IU_data, color='k', marker='o', s=5, label='TC 1 measurements')  # 's' controls the marker size
    plt.scatter(TC2_data, IU_data, color='k', marker='^', s=5, label='TC 2 measurements')
    plt.xlabel("$T$ [$^{\\circ}$C]")
    plt.ylabel("Camera counts [IU]")
    plt.legend()
    plt.savefig("/home/jpe/VKI/ELAB/2023-2024/Plots/IR_calib.pdf")
    plt.close()
    
    
    T1 = B1_optimized / np.log((R1_optimized/IU) + F1_optimized)
    T2 = B2_optimized / np.log((R2_optimized/IU) + F2_optimized)
    T_avg = ((B1_optimized / np.log((R1_optimized/IU) + F1_optimized)) + (B2_optimized / np.log((R2_optimized/IU) + F2_optimized)))/2
    
    return T_avg,T2,T1


