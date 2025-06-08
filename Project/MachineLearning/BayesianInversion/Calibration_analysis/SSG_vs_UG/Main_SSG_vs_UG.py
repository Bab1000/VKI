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

global_output_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Calibration_analysis/SSG_vs_UG"

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
#    filename_base="gamma_posteriors"
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

MAP_comparison_plot_individual(
    sm_model=sm_q,
    Pdyn_test=Pdyn_test,
    Pstat_test=Pstat_test,
    T_test=T_test,
    HF_exp_list=HF_test,
    HF_exp_err_list=0.15 * HF_test,
    global_output_path=global_output_path
)

# =====================================================================================================
 
# ----------------------
# | ERRORS COMPUTATION |
# ----------------------

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


#y_pred_ssg, y_pred_ug, error_summary = evaluate_model_comparison(
#    sm_model=sm_q,
#    Pdyn_list=Pdyn_test,
#    Pstat_list=Pstat_test,
#    T_list=T_test,
#    HF_exp_list=HF_test,
#    global_output_path=global_output_path
#)

















    #MAP_plots(
    #    sm_model=sm_q,
    #    samples=samples_PBC_SSG,
    #    targ_p=targ_p,
    #    MF_test=MF_test,
    #    Pdyn_test=Pdyn_test,
    #    T_test=T_test,
    #    T_exp_list=T_test,
    #    HF_exp_list=HF_test,
    #    HF_exp_err_list=np.array(0.15*HF_test),
    #    output_path= global_output_path + f"/MAP_predictions/PBC_SSG_{targ_p/100:.0f}mbar_predictions.jpeg",
    #    global_output_path=global_output_path
    #)
#
    #MAP_plots(
    #    sm_model=sm_q,
    #    samples=samples_PBC_UG,
    #    targ_p=targ_p,
    #    MF_test=MF_test,
    #    Pdyn_test=Pdyn_test,
    #    T_test=T_test,
    #    T_exp_list=T_test,
    #    HF_exp_list=HF_test,
    #    HF_exp_err_list=np.array(0.15*HF_test),
    #    output_path= global_output_path + f"/MAP_predictions/PBC_UG_{targ_p/100:.0f}mbar_predictions.jpeg",
    #    global_output_path=global_output_path
    #)
#
    #GN_GO_plotting(samples_PBC_SSG,targ_p,global_output_path)
#
    #GN_GO_plotting(samples_PBC_UG,targ_p,global_output_path)





    