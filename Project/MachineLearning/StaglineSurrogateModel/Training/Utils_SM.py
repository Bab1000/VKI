
import time
import sys
import threading
import numpy as np
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm
from colorama import Fore, Style, init

# -------------------------
# | DATA LOADING FUNCTION |
# -------------------------
def load_data(y_file, x_file, X_list, Y_list, desc="Loading Data"):
    """Loads input and output data from files while showing a tqdm progress bar."""
    try:
        with open(y_file) as y_f, open(x_file) as x_f:
            y_lines = y_f.readlines()
            x_lines = x_f.readlines()

            for line_index, line in tqdm(enumerate(y_lines), total=len(y_lines), desc=desc):
                if line_index > 0:
                    splitted = line.split()
                    
                    # Handling non-converged cases
                    if splitted[2] == "NotConv":
                        splitted[2] = 999
                    
                    conv = float(splitted[2])
                    
                    # If convergence is acceptable, add the data
                    if conv < 1e-2:
                        Y_list.append(float(splitted[1]))
                        x_values = x_lines[line_index].split()
                        X_list.append([float(x_values[1]), float(x_values[2]), float(x_values[3]), np.log10(float(x_values[4])), np.log10(float(x_values[5]))])

        print(Fore.GREEN + "---> [SUCCESS] Data loaded successfully!")

    except FileNotFoundError as e:
        print(Fore.RED + f"---> [ERROR] File not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Failed to load validation data: {e}")
        sys.exit(1)

# --------------------------------------
# | FUNCTION: Train Model with Logging |
# --------------------------------------
def train_model(sm_q):
    """Trains the Kriging model and prints execution time."""
    start_time = time.time()
    
    print(Fore.WHITE + "[INFO] Model training started...")
    
    sm_q.train()

    end_time = time.time()
    print(Fore.GREEN + f"[SUCCESS] Model training completed in {end_time - start_time:.2f} seconds")
