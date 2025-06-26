import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az # type: ignore
from colorama import Fore, Style, init
from Utils_Main_SSG_vs_UG import *

# Initialize colorama for colored terminal output
init(autoreset=True)

# ========================================================================
 
# ------------------------------------
# | STAGLINE SURROGATE MODEL LOADING |
# ------------------------------------

model_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/KRGModel/SM_Gamma_Log"

sm_q = LoadModel(model_path)

# ========================================================================
 
# -----------------------------
# | EXPERIMENTAL DATA LOADING |
# -----------------------------

# Path to test campaign data
CSV_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Data/tests-mass-flow-Justin.xlsx"

# Columns to check for wrong data
columns_to_check = ["Pressure[mbar]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)", "HeatFlux(HS50mm)[kW/m2]"]

Pstat,massflow,power,HF,off_set_HF,Pdyn,off_set_Pdyn,temperature = CSVReader(CSV_path,columns_to_check)

global_output_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Calibration_analysis/Models_comparison"

# ========================================================================
 
# ----------------------------------
# | GAMMA DISTRIBUTIONS PLOTTING   |
# ----------------------------------

print(Fore.BLUE + "[STEP] Computing gammas distributions")

idata_dict = {
    1500: (az.from_netcdf("../Gamma_models/PBC_SSG_15mbar.nc"), az.from_netcdf("../Gamma_models/PBC_UG_15mbar.nc")),
    5000: (az.from_netcdf("../Gamma_models/PBC_SSG_50mbar.nc"), az.from_netcdf("../Gamma_models/PBC_UG_50mbar.nc")),
    10000: (az.from_netcdf("../Gamma_models/PBC_SSG_100mbar.nc"), az.from_netcdf("../Gamma_models/PBC_UG_100mbar.nc"))
}

#plot_posteriors_individual(
#    idata_dict=idata_dict,
#    pressure_labels=[1500, 5000, 10000],
#    output_path=global_output_path,
#    filename_base="Full_gamma_posteriors"
#)


# =====================================================================================================
 
# -------------------
# | MAP PLOTTING UG |
# -------------------

print(Fore.BLUE + "[STEP] Computing Heat Fluxes")

# Target static pressure 
# ----------------------
target_pressure = [1500,5000,10000]
targ_mf = 16
tolerance_press = 300
tolerance_mf = 3

# Gathering data for bayesian inversion
# -------------------------------------

Pdyn_test = []
Pstat_test = []
T_test = []
HF_test = []
off_set_Pdyn_test = []
off_set_HF_test = []
MF_test = []

for targ_p in target_pressure:

    for i,pstat in enumerate(Pstat):

        if pstat >= targ_p - tolerance_press and pstat <= targ_p + tolerance_press:

            if massflow[i] >= targ_mf - tolerance_mf and massflow[i] <= targ_mf + tolerance_mf:

                Pdyn_test.append(Pdyn[i])
                Pstat_test.append(Pstat[i])
                T_test.append(temperature[i])
                HF_test.append(HF[i])
                off_set_Pdyn_test.append(off_set_Pdyn[i])
                off_set_HF_test.append(off_set_HF[i])
                MF_test.append(massflow[i])

HF_test = np.array(HF_test)/1000

#MAP_comparison_plot_individual(
#    sm_model=sm_q,
#    Pdyn_test=Pdyn_test,
#    Pstat_test=Pstat_test,
#    T_test=T_test,
#    HF_exp_list=HF_test,
#    HF_exp_err_list=0.15 * HF_test,
#    global_output_path=global_output_path
#)

# =====================================================================================================
 
# ----------------------
# | ERRORS COMPUTATION |
# ----------------------

# Only for panerai and viladegut models since the errors for PBC_SSG and PBC_UG were computed in the folder SSG_vs_UG

# MAPE PCB_SSG : 8.48%
# MAPE PCB_SSG : 10.27%
# Panerai → MAPE = 27.32 %, MAE = 538.965 kW/m²
# Viladegut → MAPE = 23.91 %, MAE = 446.966 kW/m²

print(Fore.BLUE + "[STEP] Errors computation")

# Target static pressure 
# ----------------------
target_pressure = [1500,5000,10000]
tolerance_press = 300
tolerance_mf = 3

# Gathering data for bayesian inversion
# -------------------------------------
Pdyn_test = []
Pstat_test = []
T_test = []
HF_test = []
off_set_Pdyn_test = []
off_set_HF_test = []
MF_test = []

for targ_p in target_pressure:

    for i,pstat in enumerate(Pstat):

        if pstat >= targ_p - tolerance_press and pstat <= targ_p + tolerance_press:

            Pdyn_test.append(Pdyn[i])
            Pstat_test.append(Pstat[i])
            T_test.append(temperature[i])
            HF_test.append(HF[i])
            off_set_Pdyn_test.append(off_set_Pdyn[i])
            off_set_HF_test.append(off_set_HF[i])
            MF_test.append(massflow[i])

HF_test = np.array(HF_test)/1000


results = evaluate_pan_and_vil_models(
    sm_model=sm_q,
    Pdyn_list=Pdyn_test,
    Pstat_list=Pstat_test,
    T_list=T_test,
    HF_exp_list=HF_test
)





# Gamma model values
pressures_new_model = [15, 50, 100]
gamma_new_model = [0.01713, 0.03265, 0.14343]

pressures_panerai = [15, 50, 100]
gamma_panerai = [0.1, 0.01, 0.005]

pressures_viladegut = [15, 50, 100]
gamma_viladegut = [0.0882, 0.02661, 0.0096]

# Plot
plt.figure(figsize=(8, 6))

plt.plot(pressures_new_model, gamma_new_model, marker='o', label='New model', color='tab:blue')
plt.plot(pressures_panerai, gamma_panerai, marker='s', label='Panerai', color='tab:orange')
plt.plot(pressures_viladegut, gamma_viladegut, marker='^', label='Viladegut & Chazot', color='tab:green')

plt.xlabel('Static Pressure $p_c$ [mbar]', fontsize=16)
plt.ylabel(r'Catalytic Efficiency $\gamma$', fontsize=16)
plt.title(r'Comparison of $\gamma$ Models vs Pressure', fontsize=16)
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save plot
plt.savefig(global_output_path + "/Model_comparison.png", dpi=300)

