import numpy as np
from colorama import Fore, Style, init
from smt.surrogate_models import KRG
import sys
import pandas as pd


# Initialize colorama for colored terminal output
init(autoreset=True)

def LoadModel(model_path):

    print(Fore.BLUE + "[STEP] Loading the surrogate model")

    print(Fore.WHITE + f"---> [INFO] Loading Model: '{model_path}' ...")
    try:
        sm_q = KRG.load(model_path)
        sm_q.options["print_global"] = False
        print(Fore.GREEN + f"---> [SUCCESS] Model '{model_path}' loaded successfully!")
        return sm_q
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Failed to load model '{model_path}': {e}")
        sys.exit(1)

def FuncBuild(sm_q,Pdyn,Pc,Tinlet,Gamma_vec):
    
    gammaN, gammaO = Gamma_vec[0], Gamma_vec[1]
    
    XV = np.array([Pdyn,Pc,Tinlet,gammaN,gammaO]).reshape(1, -1)

    try:
        # Y predictions
        YV_K= sm_q.predict_values(XV)/1000
        # Get prediction uncertainty (variance)
        #YV_K_var = sm_q.predict_variances(XV) / 1000000  # Convert from W² to kW²
        #YV_K_std = np.sqrt(YV_K_var)  # Convert variance to standard deviation
        return YV_K
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Prediction failed: {e}")
        sys.exit(1)

def FuncBuild_Full(sm_q,Pdyn,Pc,Tinlet,Gamma_vec):
    
    gammaN, gammaO = np.ones((len(Pdyn),))*np.log10(Gamma_vec[0]), np.ones((len(Pdyn),))*np.log10(Gamma_vec[1])
    
    XV = np.array([Pdyn,Pc,Tinlet,gammaN,gammaO]).reshape(len(Pdyn), -1)

    try:
        # Y predictions
        YV_K= sm_q.predict_values(XV)/1000
        # Get prediction uncertainty (variance)
        #YV_K_var = sm_q.predict_variances(XV) / 1000000  # Convert from W² to kW²
        #YV_K_std = np.sqrt(YV_K_var)  # Convert variance to standard deviation
        return YV_K
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Prediction failed: {e}")
        sys.exit(1)

def CSVReader(CSV_path,columns_to_check):
    # Reading the CSV file with all the data
    df = pd.read_excel(CSV_path, engine="openpyxl")

    # Replace "NA" strings with actual NaN values
    df.replace("NA", np.nan, inplace=True)

    # Drop rows where any of the specified columns contain NaN
    df_cleaned = df.dropna(subset=columns_to_check)

    # Gathering the cleaned data
    pressure = df_cleaned["Pressure[mbar]"].dropna().tolist()
    massflow = df_cleaned["massflow [g/s]"].dropna().tolist()
    power = df_cleaned["Power[kW]"].dropna().tolist()
    heat_flux = df_cleaned["HeatFlux(HS50mm)[kW/m2]"].dropna().tolist()
    pitot = df_cleaned["Pitot[Pa]"].dropna().tolist()
    temperature = df_cleaned["T [K] (x = 375mm, r = 0mm)"].dropna().tolist()

    return pressure,massflow,power,heat_flux,pitot,temperature

def ExpDataLoading(CSV_path,press):

    # Define the columns to check for NaN values
    columns_to_check = ["Pressure[mbar]", "massflow [g/s]", "Power[kW]", "HeatFlux(HS50mm)[kW/m2]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)"]

    # Load data from CSV
    pressure, massflow, power, heat_flux, pitot, temperature = CSVReader(CSV_path, columns_to_check)

    # Define target mdot values
    target_mdot_values = [10, 16, 20]

    # Define tolerance range for filtering
    range_tolerance = 3

    # Dictionary to store experimental heat flux and temperature for each mdot
    exp_HF_data = {mdot: [] for mdot in target_mdot_values}
    T_mdot_data = {mdot: [] for mdot in target_mdot_values}
    pitot_mdot_data = {mdot: [] for mdot in target_mdot_values}

    # Filter experimental data based on target pressure and mdot values
    for i, exp_press in enumerate(pressure):
        if press - range_tolerance <= exp_press <= press + range_tolerance:
            for mdot in target_mdot_values:
                if mdot - range_tolerance <= massflow[i] <= mdot + range_tolerance:
                    exp_HF_data[mdot].append(heat_flux[i])
                    T_mdot_data[mdot].append(temperature[i])
                    pitot_mdot_data[mdot].append(pitot[i])

    # Convert lists to NumPy arrays for better handling
    for mdot in target_mdot_values:
        exp_HF_data[mdot] = np.array(exp_HF_data[mdot])
        T_mdot_data[mdot] = np.array(T_mdot_data[mdot])
        pitot_mdot_data[mdot] = np.array(pitot_mdot_data[mdot])

    # Sorting based on exp_HF_data while keeping T_mdot_data aligned
    for mdot in exp_HF_data:
        sorted_indices = np.argsort(exp_HF_data[mdot])  # Get sorted indices
        exp_HF_data[mdot] = exp_HF_data[mdot][sorted_indices]  # Sort exp_HF_data
        T_mdot_data[mdot] = T_mdot_data[mdot][sorted_indices]  # Align T_mdot_data
        pitot_mdot_data[mdot] = pitot_mdot_data[mdot][sorted_indices]  # Align T_mdot_data

    return exp_HF_data, T_mdot_data, pitot_mdot_data, target_mdot_values

def FullDataLoading(CSV_path):
    # Define the columns to check for NaN values
    columns_to_check = ["Pressure[mbar]", "massflow [g/s]", "Power[kW]", "HeatFlux(HS50mm)[kW/m2]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)"]

    # Load data from CSV
    pressure, massflow, power, heat_flux, pitot, temperature = CSVReader(CSV_path, columns_to_check)

    return pressure, heat_flux, pitot, temperature

