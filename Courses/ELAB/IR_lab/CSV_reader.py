import os
import csv
import pandas as pd
import pdb

def CSV_reader(dic,folder):
    # Iterate over all files in the specified directory
    
    for filename in os.listdir(folder):
        # Check if the file is a CSV
        if filename.endswith('.csv'):
            # Create the full file path
            file_path = os.path.join(folder, filename)
            
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path, sep=';', header=None).to_numpy()      # Change the separator if necessary
            
            # Store the DataFrame in the dictionary with the filename (without extension) as the key
            key_name = os.path.splitext(filename)[0]  # Get the filename without extension
            dic[key_name] = df


    # Step 1: Get the list of keys sorted alphabetically
    sorted_keys = sorted(dic.keys())

    # Step 2: Create an empty dictionary for storing the sorted data
    sorted_dic = {}

    # Step 3: Iterate over the sorted keys and populate the new dictionary
    for key in sorted_keys:
        sorted_dic[key] = dic[key]

    return sorted_dic