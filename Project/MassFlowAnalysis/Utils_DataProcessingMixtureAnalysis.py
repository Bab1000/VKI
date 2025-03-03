import re
import numpy as np
import pandas as pd



# Function to extract mdot from the file path
def extract_mdot(filepath):
    match = re.search(r"mdot=([\d\.]+)", filepath)
    return float(match.group(1)) if match else None

# Function to extract temperature T from the file path
def extract_temperature(filepath):
    match = re.search(r"T=([\d\.]+)", filepath)
    return float(match.group(1)) if match else None

# Function to extract q_tot from the file content
def extract_q_tot(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'q_tot \[W/m2\]:\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                return float(match.group(1))
    return None

# Function to round mdot values to the closest valid option (10, 16, or 20)
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