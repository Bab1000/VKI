import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
from colorama import Fore, Style, init
from smt.surrogate_models import KRG
from Utils_SM import *

# Initialize colorama for colored terminal output
init(autoreset=True)

# ===========================================================================================

# --------------------------
# | LOAD PRE-TRAINED MODEL |
# --------------------------
print(Fore.BLUE + "[STEP] Loading the surrogate model")

model_folder_name = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Models"

model_path = model_folder_name + "/SM_Gamma_Log"

print(Fore.WHITE + f"---> [INFO] Loading Model: '{model_path}' ...")
try:
    sm_q = KRG.load(model_path)
    print(Fore.GREEN + f"---> [SUCCESS] Model '{model_path}' loaded successfully!")
except Exception as e:
    print(Fore.RED + f"---> [ERROR] Failed to load model '{model_path}': {e}")
    sys.exit(1)

# ===========================================================================================

# ------------------------------
# | INPUTS FOR SURROGATE MODEL |
# ------------------------------
print(Fore.BLUE + "[STEP] Preparing data to predict")

# Data containers

temperature = np.arange(5500,7510,10)

pdyn_var = np.arange(50,350,50)

pdyn = np.full(len(temperature), 100)

pstat = np.full(len(temperature), 1500)

pstat_var = [1500, 5000, 10000, 20000]

gammaN = np.full(len(temperature), -1)

gamma0 = np.full(len(temperature), -1)

# ===========================================================================================

# ----------------------------------
# | MODEL PREDICTION pdyn VARIATION|
# ----------------------------------


print(Fore.BLUE + "[STEP] Model predictions for pdyn variations")
print(Fore.WHITE + "---> [INFO] Running predictions on data ...")
YV_pdyn_list = []
for pdyn_i in pdyn_var:
    try:
        # Create list of inputs
        pdyn_full = np.full(len(temperature), pdyn_i)
        XV = np.column_stack((pdyn_full,pstat,temperature,gammaN,gamma0))
        # Convert lists to NumPy arrays
        XV = np.array(XV)
        # Y predictions
        YV_K= sm_q.predict_values(XV)/1000
        YV_pdyn_list.append(YV_K)
        # Get prediction uncertainty (variance)
        YV_K_var = sm_q.predict_variances(XV) / 1000000  # Convert from W² to kW²
        YV_K_std = np.sqrt(YV_K_var)  # Convert variance to standard deviation
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Prediction failed: {e}")
        sys.exit(1)
    
print(Fore.GREEN + "---> [SUCCESS] Model Predictions successful for pdyn variations!")

# ===========================================================================================

# ----------------------------------
# | MODEL PREDICTION pstat VARIATION|
# ----------------------------------


print(Fore.BLUE + "[STEP] Model predictions for pstat variations")
print(Fore.WHITE + "---> [INFO] Running predictions on data ...")
YV_pstat_list = []
for pstat_i in pstat_var:
    try:
        # Create list of inputs
        pstat_full = np.full(len(temperature), pstat_i)
        XV = np.column_stack((pdyn,pstat_full,temperature,gammaN,gamma0))
        # Convert lists to NumPy arrays
        XV = np.array(XV)
        # Y predictions
        YV_K= sm_q.predict_values(XV)/1000
        YV_pstat_list.append(YV_K)
        # Get prediction uncertainty (variance)
        YV_K_var = sm_q.predict_variances(XV) / 1000000  # Convert from W² to kW²
        YV_K_std = np.sqrt(YV_K_var)  # Convert variance to standard deviation
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Prediction failed: {e}")
        sys.exit(1)
    
print(Fore.GREEN + "---> [SUCCESS] Model Predictions successful for pstat variations!")

# ===========================================================================================

# ----------------
# | PLOT RESULTS |
# ----------------

print(Fore.BLUE + "[STEP] Plotting the results")

res_folder_name = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Training/Results"

print(Fore.WHITE + "---> [INFO] Generating visualization ...")

# Apply a clean style
plt.style.use("seaborn-v0_8-whitegrid")
mpl.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

# Create the figure
plt.figure(figsize=(7, 6))

# Plot each line for varying pdyn values
for i in range(len(YV_pdyn_list)):
    label = f"Pdyn: {pdyn_var[i]} Pa"
    plt.plot(temperature, YV_pdyn_list[i], label=label, linewidth=2)

# Add vertical line at T = 6000 K
plt.axvline(x=6000, color='red', linestyle='--', linewidth=2)
#plt.axvline(x=7000, color='red', linestyle='--', linewidth=2)
# Get axis limits
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
# Compute X positions (25% and 75%)
x_air5 = xmin + 0.25 * (xmax - xmin)
x_air11 = xmin + 0.75 * (xmax - xmin)
# Y position: slightly above the bottom (to avoid overlap with x-axis)
y_text = ymin + 0.02 * (ymax - ymin)  # 2% above the bottom
# Add the labels
plt.text(5750, y_text, "air_5", fontsize=13, fontweight='bold', color='black', ha='center', va='bottom')
plt.text(6750, y_text, "air_11", fontsize=13, fontweight='bold', color='black', ha='center', va='bottom')
# Axis labels and title
plt.xlabel(r"$\it{Temperature}$ [K]")
plt.ylabel(r"$q_{wall}$ predicted [kW]")
# Display the legend with background and shadow
plt.legend(loc="upper left")
# Ensure layout fits well
plt.tight_layout()
# Save the figure in high resolution
plt.savefig(res_folder_name + "/MixtureValidation_pdynVariations.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
plt.close()

# Create the plot
plt.figure(figsize=(7, 6))
# Plot each line for varying pstat values
for i in range(len(YV_pstat_list)):
    label = f"Pc: {pstat_var[i]/100:.0f} mbar"
    plt.plot(temperature, YV_pstat_list[i], label=label, linewidth=2)

# Add vertical line at T = 6000 K
plt.axvline(x=6000, color='red', linestyle='--', linewidth=2)
#plt.axvline(x=7000, color='red', linestyle='--', linewidth=2)
# Get axis limits
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
# Compute X positions (25% and 75%)
x_air5 = xmin + 0.25 * (xmax - xmin)
x_air11 = xmin + 0.75 * (xmax - xmin)
# Y position: slightly above the bottom (to avoid overlap with x-axis)
y_text = ymin + 0.02 * (ymax - ymin)  # 2% above the bottom
# Add the labels
plt.text(5750, y_text, "air_5", fontsize=13, fontweight='bold', color='black', ha='center', va='bottom')
plt.text(6750, y_text, "air_11", fontsize=13, fontweight='bold', color='black', ha='center', va='bottom')
# Axis labels and title
plt.xlabel(r"$\it{Temperature}$ [K]")
plt.ylabel(r"$q_{wall}$ predicted [kW]")
# Display the legend with background and shadow
plt.legend(loc="upper left")
# Ensure layout fits well
plt.tight_layout()
# Save the figure in high resolution
plt.savefig(res_folder_name + "/MixtureValidation_pstatVariations.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
plt.close()

print(Fore.GREEN + f"---> [SUCCESS] Results successfully ploted in {res_folder_name + "/ModelPerformance.pdf"}")
