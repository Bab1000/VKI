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

CSV_name = f"/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Res_tuning_Pstat_Stagline_Single_dx.csv"

df = pd.read_csv(CSV_name)

index = set(df["index"])

# Créer l'ensemble des indices valides (comme tu filtres XT)
valid_idx = [idx for idx, line in enumerate(XT) if idx in index and line[2] > 7000]

# Filtrer df en conséquence
df_filtered = df[df["index"].isin(valid_idx)].reset_index(drop=True)

# Extraire les colonnes filtrées
T = df_filtered["Tinlet"].values
Pstat = df_filtered["Pstat"].values / 100
dx_diff = df_filtered["dx_diff"].values
dx_conv = df_filtered["dx_conv"].values

inputs = []
for idx, line in enumerate(XT): 
    if line[2] > 7000 :
        inputs.append(line)

inputs = [line for idx, line in enumerate(XT) if idx in index]

gamma_n = []
gamma_o = []
for idx,input in enumerate(inputs):
    gamma_n.append(input[3])
    gamma_o.append(input[4])

gamma_n = np.array(gamma_n)
gamma_o = np.array(gamma_o)

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
    Y=dx_diff
)

results_dx_conv = train_compare_models(
    T=T,
    Pstat=Pstat,
    Y=dx_conv
)

# ======================================================================================================================

# Plotting

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(T, Pstat, np.log10(dx_diff))

ax.set_xlabel('T [K]')
ax.set_ylabel('Pstat [mbar]')
ax.set_zlabel('log10(dx_diff)')

plt.tight_layout()
plt.savefig(f"/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Plots/Res_dx_diff.jpeg",format='jpeg', dpi=300, bbox_inches='tight')
#plt.show()
plt.close()


fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(T, Pstat, np.log10(dx_conv))

ax.set_xlabel('T [K]')
ax.set_ylabel('Pstat [mbar]')
ax.set_zlabel('log10(dx_conv)')

plt.tight_layout()
plt.savefig(f"/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Plots/Res_dx_conv.jpeg",format='jpeg', dpi=300, bbox_inches='tight')
#plt.show()
plt.close()


