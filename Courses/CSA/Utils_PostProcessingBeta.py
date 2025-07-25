import re
import numpy as np
import pandas as pd

def extract_mdot(filepath):
    match = re.search(r"mdot=([\d\.]+)", filepath)
    return float(match.group(1)) if match else None

def extract_temperature(filepath):
    match = re.search(r"T=([\d\.]+)", filepath)
    return float(match.group(1)) if match else None

def extract_pressure(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r"Pressure\s+\[Pa\]:\s+([0-9.Ee+-]+)", line)
            if match:
                return float(match.group(1))
    return None

def extract_h_wall(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r"h_wall\s+\[J/kg\]:\s+([0-9.Ee+-]+)", line)
            if match:
                return float(match.group(1))
    return None

def extract_q_tot(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'q_tot \[W/m2\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))
    return None

def extract_h_edge(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'h_edge \[J/kg\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))
    return None

def extract_BL_edge(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'BL_thickness \[m\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))
    return None

def extract_u_edge(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'u_edge \[m/s\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))
    return None

def extract_density_edge(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'Densitiy_edge \[kg/m3\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))
    return None

def extract_R(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        if len(lines) >= 2:
            line = lines[1] 
            try:
                return float(line.strip())
            except ValueError:
                return None 
        return None
            
    

def extract_heat_transfer_coeff(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'heat_transfer_coeff \[-\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))

def round_mdot(mdot):
    return min([10, 16, 20], key=lambda x: abs(x - mdot))

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

def FlowfieldReader(flowfield_path):
    # Reading the CSV file with all the data
    df = pd.read_csv(flowfield_path, sep=r'\s+', header=None, skiprows=4)

    # Gathering the data
    mesh = df.iloc[:,0].tolist()
    U = df.iloc[:,1].tolist()
    V = df.iloc[:,2].tolist()
    T = df.iloc[:,3].tolist()
    pc = df.iloc[:,4].tolist()

    mesh = np.array(mesh)
    U = np.array(U)
    V = np.array(V)
    T = np.array(T)
    pc = np.array(pc)
    

    return mesh,U,V,T,pc