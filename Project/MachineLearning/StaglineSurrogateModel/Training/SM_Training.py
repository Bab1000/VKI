
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from colorama import Fore, Style, init
from smt.surrogate_models import KRG

from Utils_SM import *

# Initialize colorama for colored terminal output
init(autoreset=True)

# ===========================================================================================

# ----------------
# | DATA LOADING |
# ----------------

# File paths for training and validation data
data_folder_name = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Data"
train_Y_file = data_folder_name + "/YTraining.dat"
train_X_file = data_folder_name + "/Training_Xdata.dat"
valid_Y_file = data_folder_name + "/YValidation.dat"
valid_X_file = data_folder_name + "/Validation_Xdata.dat"


# Data containers
XT, YT, XV, YV = [], [], [], []

# Load training data
load_data(train_Y_file, train_X_file, XT, YT, desc="Loading Training Data")

print(Fore.GREEN + "[SUCCESS] Training data successfully loaded!")

# Load validation data
load_data(valid_Y_file, valid_X_file, XV, YV, desc="Loading Validation Data")

# Convert lists to NumPy arrays for efficient numerical operations
XT, YT, XV, YV = map(np.array, [XT, YT, XV, YV])

print(Fore.GREEN + "[SUCCESS] Validation data successfully loaded!")

# ===========================================================================================

# -------------------------------
# | MODEL DEFINITION & TRAINING |
# -------------------------------

# Initialize the Kriging model
print(Fore.WHITE + "[INFO] Initialising the surrogate model ...")
sm_q = KRG(theta0=[1e-1], corr="matern52", print_global=True)

# Set training values
print(Fore.WHITE + "[INFO] Initialising training set for the SM ...")
sm_q.set_training_values(XT, YT)

# Train the model
train_model(sm_q)

# Save trained model
model_folder_name = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Models"
model_name = "SM_Gamma_Log"
sm_q.save(model_folder_name + f"/{model_name}")
print(Fore.GREEN + f"[SUCCESS] Model saved successfully as {model_name}")

# ===========================================================================================