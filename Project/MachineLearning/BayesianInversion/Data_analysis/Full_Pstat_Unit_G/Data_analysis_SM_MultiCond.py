import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az # type: ignore
from colorama import Fore, Style, init
from Utils_Data_analysis_MCT import *

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

# ========================================================================
 
# -----------------------------
# | PROCESSING OF THE RESULTS |
# -----------------------------

global_output_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Data_analysis/Full_Pstat_Unit_G"

print(Fore.BLUE + "[STEP] Computing Heat Fluxes")

# Target static pressure 
# ----------------------
target_pressure = [1500,5000,10000]
tolerance_press = 300
tolerance_mf = 3

for targ_p in target_pressure:

    # Gathering data for bayesian inversion
    # -------------------------------------
    Pdyn_test = []
    Pstat_test = []
    T_test = []
    HF_test = []
    off_set_Pdyn_test = []
    off_set_HF_test = []
    MF_test = []



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

    # --------------------------------------------
    # | SURROGATE MODEL PREDICTIONS & PLOTTING   |
    # --------------------------------------------

    # === Getting the samples ===
    trace_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Simulations/Full_Pstat_Unit_G/UniformG_MultiCond_Pstat={targ_p/100}/Posteriors/Thinned_trace.nc"  
    idata = az.from_netcdf(trace_path)

    # === Extract posterior samples for all input variables ===
    posterior = idata.posterior
    samples = posterior.stack(sample=("chain", "draw"))

    #MAP_plots_Unit_G(
    #    sm_model=sm_q,
    #    samples=samples,
    #    targ_p=targ_p,
    #    MF_test=MF_test,
    #    Pdyn_test=Pdyn_test,
    #    T_test=T_test,
    #    T_exp_list=T_test,
    #    HF_exp_list=HF_test,
    #    HF_exp_err_list=np.array(0.15*HF_test),
    #    output_path= global_output_path + f"/MAP_HF_values/FullCond_MAP_Pstat={targ_p/100}.jpeg",
    #    global_output_path=global_output_path
    #)
#
    #G_plotting(samples,targ_p,global_output_path)
    #
    #plot_correlation(idata, targ_p, kind='scatter', show_plots=False, save_dir= global_output_path + "/Correlation")
    #
    #plot_T_analysis(samples, T_test, 0.02, global_output_path,targ_p, n_cases=10)

    plot_relative_errors_by_massflow(sm_q, samples, targ_p, T_test, HF_test, MF_test, global_output_path)



    