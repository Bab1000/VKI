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
 
Pdyn = np.linspace(40, 200, 5)               
Pstat = np.full(5, 1500)                      
T = np.linspace(6000, 8000, 5)               
GN = np.full(5, -2)                      
GO = np.full(5, -3)                        

X_vec = np.stack((Pdyn, Pstat, T, GN, GO), axis=1)

HF = sm_q.predict_values(X_vec)[:, 0]/1000  # shape (5,)

HF_error = np.array(0.15*HF)  # conversion en kW/m²

# ========================================================================
 
# -----------------------------
# | PROCESSING OF THE RESULTS |
# -----------------------------

global_output_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case_15mbar/Results_analysis"

print(Fore.BLUE + "[STEP] Computing Heat Fluxes")

# --------------------------------------------
# | SURROGATE MODEL PREDICTIONS & PLOTTING   |
# --------------------------------------------

# === Getting the samples ===
trace_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case_15mbar/Posteriors/Thinned_trace.nc"  
idata = az.from_netcdf(trace_path)

# === Extract posterior samples for all input variables ===
posterior = idata.posterior
samples = posterior.stack(sample=("chain", "draw"))

MAP_plots(
    sm_model=sm_q,
    samples=samples,
    targ_p=Pstat[0],
    Pdyn_test=Pdyn,
    T_test=T,
    T_exp_list=T,
    HF_exp_list=HF,
    HF_exp_err_list=np.array(0.15*HF),
    output_path= global_output_path + f"/MAP_HF_values/FullCond_MAP_Pstat={Pstat[0]/100}.jpeg",
    global_output_path=global_output_path
)

GN_GO_plotting(samples,Pstat[0],global_output_path)

GN_GO_separate_plotting(samples,Pstat[0],global_output_path)





