import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

from Utils_0DReactor import *

# ===========================================================================================

# ----------------
# | DATA LOADING |
# ----------------

# File paths for training and validation data
data_folder_name = "/home/jpe/VKI/Courses/PCMAF/Data"
train_Y_file = data_folder_name + "/YTraining.dat"
train_X_file = data_folder_name + "/Training_Xdata.dat"
valid_Y_file = data_folder_name + "/YValidation.dat"
valid_X_file = data_folder_name + "/Validation_Xdata.dat"


# Data containers
XT, YT, XV, YV = [], [], [], []

# Load training data
load_data(train_Y_file, train_X_file, XT, YT, desc="Loading Training Data")

print(Fore.GREEN + "[SUCCESS] Training data successfully loaded!")

CSV_name = f"/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Res_tuning_Pstat_Stagline_dx_conv_only.csv"

df = pd.read_csv(CSV_name)

index = set(df["index"])

# Extraire les colonnes filtrées
T = df["Tinlet"].values
Pstat = df["Pstat"].values / 100
dx_diff = df["dx_diff"].values
dx_conv = df["dx_conv"].values

# Calcul du masque : log10(dx) < -3 équivaut à dx < 1e-3
#mask = (np.log10(dx_diff) > -3) & (np.log10(dx_conv) > -3)

# Application du filtre
dx_diff_filtered = dx_diff
dx_conv_filtered = dx_conv
T = T
Pstat = Pstat

print(np.shape(T))

# ======================================================================================================================

# Model training

def train_compare_models(T, Pstat, Y, test_size=0.2, seed=42):
    """
    Train and compare Random Forest, SVR, and GPR on (T, Pstat, gamma_N, gamma_O) → log10(Y)
    """

    # Prepare data
    X = np.column_stack((T, Pstat))
    y = Y

    # Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)

    results = {}

    # === Random Forest ===
    rf_pipeline = make_pipeline(
        StandardScaler(),
        RandomForestRegressor(n_estimators=300, random_state=seed)
    )

    rf_pipeline.fit(X_train, y_train)
    y_pred_rf = rf_pipeline.predict(X_test)

    mae = mean_absolute_error(y_test, y_pred_rf)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred_rf))
    r2 = r2_score(y_test, y_pred_rf)
    accuracy_like = 1 - (mae / np.mean(y_test))

    results['Random Forest'] = {
        'MAE': mae,
        'RMSE': rmse,
        'R2': r2,
        'Accuracy': accuracy_like
    }

    # === SVR ===
    svr = make_pipeline(
        StandardScaler(),
        SVR(kernel='rbf', C=10, epsilon=0.01)
    )
    svr.fit(X_train, y_train)
    y_pred_svr = svr.predict(X_test)

    mae = mean_absolute_error(y_test, y_pred_svr)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred_svr))
    r2 = r2_score(y_test, y_pred_svr)
    accuracy_like = 1 - (mae / np.mean(y_test))

    results['SVR'] = {
        'MAE': mae,
        'RMSE': rmse,
        'R2': r2,
        'Accuracy': accuracy_like
    }

    # === Gaussian Process (Kriging) ===
    kernel = C(1.0, (1e-2, 1e3)) * RBF(length_scale=1.0)
    gpr = GaussianProcessRegressor(kernel=kernel, alpha=1e-6, normalize_y=True, random_state=seed)
    gpr.fit(X_train, y_train)
    y_pred_gpr = gpr.predict(X_test)

    mae = mean_absolute_error(y_test, y_pred_gpr)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred_gpr))
    r2 = r2_score(y_test, y_pred_gpr)
    accuracy_like = 1 - (mae / np.mean(y_test))

    results['GPR'] = {
        'MAE': mae,
        'RMSE': rmse,
        'R2': r2,
        'Accuracy': accuracy_like
    }

    # === Display results ===
    print("\n🔍 Comparison of models on test set (log10(Y)):\n")
    for name, res in results.items():
        print(f"{name}:")
        print(f"  → MAE  = {res['MAE']:.4f}")
        print(f"  → RMSE = {res['RMSE']:.4f}")
        print(f"  → R²   = {res['R2']:.4f}")
        print(f"  → Accuracy   = {res['Accuracy']:.4f}\n")

    return results

results_dx_diff = train_compare_models(
    T=T,
    Pstat=Pstat,
    Y=dx_diff_filtered
)

results_dx_conv = train_compare_models(
    T=T,
    Pstat=Pstat,
    Y=dx_conv_filtered
)

# ======================================================================================================================
# Plotting avec colorbar

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

log_dx_diff = np.log10(dx_diff_filtered)
sc1 = ax.scatter(T, Pstat, log_dx_diff, c=log_dx_diff, cmap='viridis')
cb1 = plt.colorbar(sc1, ax=ax, pad=0.1)
cb1.set_label('log10(dx_diff)', fontsize=12)

ax.set_xlabel('T [K]')
ax.set_ylabel('Pstat [mbar]')
ax.set_zlabel('log10(dx_diff)')

plt.tight_layout()
plt.savefig("/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Plots/Res_dx_diff_colored.jpeg",
            format='jpeg', dpi=300, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

log_dx_conv = np.log10(dx_conv_filtered)
sc2 = ax.scatter(T, Pstat, log_dx_conv, c=log_dx_conv, cmap='plasma')
cb2 = plt.colorbar(sc2, ax=ax, pad=0.1)
cb2.set_label('log10(dx_conv)', fontsize=12)

ax.set_xlabel('T [K]')
ax.set_ylabel('Pstat [mbar]')
ax.set_zlabel('log10(dx_conv)')

plt.tight_layout()
plt.savefig("/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Plots/Res_dx_conv_colored.jpeg",
            format='jpeg', dpi=300, bbox_inches='tight')
plt.close()


import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# === Données ===
X_plot = np.column_stack((T, Pstat))  # shape (N, 2)
Z_plot = np.log10(dx_conv_filtered)   # shape (N,)

# === Création d’une grille régulière ===
print(f"log10(dx_diff) min = {Z_plot.min():.4f}, max = {Z_plot.max():.4f}")
T_lin = np.linspace(T.min(), T.max(), 200)
Pstat_lin = np.linspace(Pstat.min(), Pstat.max(), 200)
T_grid, Pstat_grid = np.meshgrid(T_lin, Pstat_lin)

# === Interpolation sur la grille ===
Z_grid = griddata(X_plot, Z_plot, (T_grid, Pstat_grid), method='linear')

# === Affichage ===
plt.figure(figsize=(10, 7))
contour = plt.contourf(T_grid, Pstat_grid, Z_grid, levels=100, cmap='viridis')
cb = plt.colorbar(contour)
cb.set_label('log10(dx_diff)', fontsize=12)

plt.xlabel('T [K]')
plt.ylabel('Pstat [mbar]')
plt.title('Colormap log10(dx_diff)')
plt.tight_layout()
plt.savefig("/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Plots/Colormap_dx.jpeg",
            format='jpeg', dpi=300, bbox_inches='tight')
plt.close()
