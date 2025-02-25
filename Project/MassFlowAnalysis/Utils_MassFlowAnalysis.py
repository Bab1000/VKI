import pandas as pd

def CSVReader(CSV_path,columns_to_check):
    # Reading the CSV file with all the data
    df = pd.read_excel(CSV_path, engine="openpyxl")

    # Drop rows where any of these columns contain NaN
    df_cleaned = df.dropna(subset=columns_to_check)

    # Gathering the data
    pressure = df["Pressure[mbar]"].tolist()
    massflow = df["massflow [g/s]"].tolist()
    power = df["Power[kW]"].tolist()
    pitot = df["Pitot[Pa]"].tolist()
    temperature = df["T [K] (x = 375mm, r = 0mm)"].tolist()

    return pressure,massflow,power,pitot,temperature