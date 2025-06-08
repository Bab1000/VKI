import sys
import numpy as np
import matplotlib.pyplot as plt
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

# ----------------
# | LOADING DATA |
# ----------------
print(Fore.BLUE + "[STEP] Loading the data")
data_folder_name = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Data"
validation_y_file = data_folder_name + "/YValidation.dat"
validation_x_file = data_folder_name + "/Validation_Xdata.dat"

# Data containers
XV, YV = [], []

# Load validation data
load_data(validation_y_file, validation_x_file, XV, YV,desc="---> [INFO] Loading validation data")

# Convert lists to NumPy arrays
XV, YV = np.array(XV), np.array(YV) / 1000  # Convert to kW
print(f"Total training points:{len(XV)}")

# ===========================================================================================

# --------------------
# | MODEL PREDICTION |
# --------------------
print(Fore.BLUE + "[STEP] Model predictions")
print(Fore.WHITE + "---> [INFO] Running predictions on validation data ...")
try:
    # Y predictions
    YV_K= sm_q.predict_values(XV)/1000
    # Get prediction uncertainty (variance)
    YV_K_var = sm_q.predict_variances(XV) / 1000000  # Convert from W² to kW²
    YV_K_std = np.sqrt(YV_K_var)  # Convert variance to standard deviation
    print(Fore.GREEN + "---> [SUCCESS] Model Predictions successful!")
except Exception as e:
    print(Fore.RED + f"---> [ERROR] Prediction failed: {e}")
    sys.exit(1)

# Compute NRMSE
err = 0
for i in range(len(YV_K)):
    err += (YV[i] - YV_K[i])**2  # Squared error accumulation

nrmse = np.sqrt(err/len(YV_K)) * 100 / (np.max(YV_K) - np.min(YV_K))


# Ensure `rmse` is a float before printing
print(Fore.YELLOW + f"---> [RESULT] NRMSE Error: {float(nrmse[0]):.4f}%")

# ===========================================================================================

# ----------------
# | PLOT RESULTS |
# ----------------

print(Fore.BLUE + "[STEP] Plotting the results")

res_folder_name = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Training/Results"

print(Fore.WHITE + "---> [INFO] Generating visualization ...")

# Ensure data is sorted by X-values for smooth confidence intervals
sorted_indices = np.argsort(YV)  # Sort based on the 4th feature (T [K])
Y_sorted = YV[sorted_indices].flatten()  # Extract and flatten feature values
YV_K_sorted = YV_K[sorted_indices].flatten()  # Flatten predictions
YV_K_std_sorted = np.sqrt(YV_K_var[sorted_indices]).flatten()  # Flatten standard deviation
# Define confidence intervals (99% corresponds to ±3 standard deviations)
YV_K_lower = YV_K_sorted - 2 * YV_K_std_sorted
YV_K_upper = YV_K_sorted + 2 * YV_K_std_sorted

# Create the plot
plt.figure(figsize=(6, 6))
# Plot the confidence interval with smooth shading
plt.plot([np.min(YV), np.max(YV)], [np.min(YV), np.max(YV)], 'r--', label="Perfect Fit")
#plt.fill_between(
#    Y_sorted, YV_K_lower, YV_K_upper,
#    color="orange", alpha=0.4, label=r"95% Confidence Interval $[kW^2]$"
#)
# Scatter plot for actual and predicted values
plt.scatter(YV, YV_K, label="Stagline data vs predictions", alpha=0.7)
# Add labels and title
plt.legend(loc="upper left")
plt.xlabel(r"$q_{wall}$ from Stagline [kW]", fontsize=14)
plt.ylabel(r"$q_{wall}$ predicted [kW]", fontsize=14)
# Save figure in high resolution
plt.savefig(res_folder_name + "/ModelPerformance.jpeg", format='jpeg', dpi=300, bbox_inches='tight')

plt.figure(figsize=(7, 6))

# Identity line
plt.plot([np.min(YV), np.max(YV)], [np.min(YV), np.max(YV)], 'r--', linewidth=2, label="Perfect Fit")

# Scatter plot
plt.scatter(YV, YV_K, label=r"Stagline vs Surrogate", alpha=0.7)

# Axis labels
plt.xlabel(r"$q_{\mathrm{wall}}$ from STAGLINE [kW]", fontsize=14)
plt.ylabel(r"$q_{\mathrm{wall}}$ predicted [kW]", fontsize=14)

# Legend
plt.legend(loc="upper left", fontsize=12, frameon=False, shadow=False)

# Grid
plt.grid(True)

# Match spines color to grid color
ax = plt.gca()
grid_color = ax.yaxis.get_gridlines()[0].get_color()
for spine in ax.spines.values():
    spine.set_edgecolor(grid_color)
    spine.set_linewidth(1.0)

# Save
plt.tight_layout()
plt.savefig(res_folder_name + "/ModelPerformance_HomogeneousStyle.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
plt.close()

print(Fore.GREEN + f"---> [SUCCESS] Results successfully plotted in {res_folder_name + '/ModelPerformance.jpeg'}")
