import numpy as np
import pandas as pd

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
    pitot = df_cleaned["Pitot[Pa]"].dropna().tolist()
    temperature = df_cleaned["T [K] (x = 375mm, r = 0mm)"].dropna().tolist()

    return pressure,massflow,power,pitot,temperature