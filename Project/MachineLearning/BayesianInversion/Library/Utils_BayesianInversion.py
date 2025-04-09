import pandas as pd
import numpy as np
from smt.surrogate_models import KRG
import sys
from colorama import Fore, Style, init

# Initialize colorama for colored terminal output
init(autoreset=True)

def CSVReader(CSV_path,columns_to_check):
    try:

        print(Fore.BLUE + "[STEP] Loading data from CSV")

        # Reading the CSV file with all the data
        df = pd.read_excel(CSV_path, engine="openpyxl")

        print(Fore.WHITE + "---> [INFO] Treating data ...")
        # Replace "NA" strings with actual NaN values
        df.replace("NA", np.nan, inplace=True)

        # Drop rows where any of the specified columns contain NaN
        df_cleaned = df.dropna(subset=columns_to_check)

        print(Fore.WHITE + "---> [INFO] Preparing the data for external use ...")
        # Gathering the cleaned data
        pressure = np.array(df_cleaned["Pressure[mbar]"].dropna()) * 100 # [mbar] --> [Pa]
        massflow = np.array(df_cleaned["massflow [g/s]"].dropna())
        power = np.array(df_cleaned["Power[kW]"].dropna())
        heat_flux = np.array(df_cleaned["HeatFlux(HS50mm)[kW/m2]"].dropna()) * 1000 # [kW/m2] --> [W/m2]
        off_set_heat_flux = np.array(df_cleaned["offsetHF[kW/m2]"].dropna()) * 1000 # [kW/m2] --> [W/m2]
        pitot = np.array(df_cleaned["Pitot[Pa]"].dropna())
        off_set_pitot = np.array(df_cleaned["offsetPitot[Pa]"].dropna())
        temperature = np.array(df_cleaned["T [K] (x = 375mm, r = 0mm)"].dropna())

        # Ensure + values in off_sets
        off_set_heat_flux = abs(off_set_heat_flux)
        off_set_pitot = abs(off_set_pitot)

        print(Fore.GREEN + "---> [SUCCESS] Data loaded successfully !")

        return pressure,massflow,power,heat_flux,off_set_heat_flux,pitot,off_set_pitot,temperature

    except Exception as e:
        print(Fore.RED + f"[ERROR] Data loading failed: {e}")

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