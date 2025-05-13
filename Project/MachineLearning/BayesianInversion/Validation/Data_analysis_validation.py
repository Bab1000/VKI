import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az # type: ignore
from colorama import Fore, Style, init

from Library.Utils_BayesianInversion import *
from Utils_Data_analysis_validation import *

# Initialize colorama for colored terminal output
init(autoreset=True)

# ========================================================================
 
# ------------------------------------
# | STAGLINE SURROGATE MODEL LOADING |
# ------------------------------------

model_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/KRGModel/SM_Gamma_Log"

sm_q = LoadModel(model_path)

# ========================================================================


# --------------------------------------------
# | SURROGATE MODEL PREDICTIONS & PLOTTING   |
# --------------------------------------------

Pdyn = np.linspace(40, 200, 5)               
Pstat = np.full(5, 1500)                      
T = np.linspace(6000, 8000, 5)               
GN = np.full(5, -2)                      
GO = np.full(5, -3)                        

X_vec = np.stack((Pdyn, Pstat, T, GN, GO), axis=1)

HF = sm_q.predict_values(X_vec)[:, 0]/1000  # shape (5,)

HF_error = np.array(0.15*HF)  # conversion en kW/m²

# === Getting the samples ===
trace_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case/Posteriors/Thinned_trace.nc"  
idata = az.from_netcdf(trace_path)

# === Extract posterior samples for all input variables ===
posterior = idata.posterior
samples = posterior.stack(sample=("chain", "draw"))

plot_joint_GN_GO(samples)

Mean_plots(
    sm_model=sm_q,
    T_test=T,
    HF_test=HF,
    HF_error=HF_error,
    Pdyn_test=Pdyn,
    Pstat_test=Pstat,
    trace_path=trace_path,
    output_path=f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case/Plots/Mean_value_validation.jpeg"
)

MAP_plots(
    sm_model=sm_q,
    samples=samples,
    Pdyn_test=Pdyn,
    T_test=T,
    T_exp_list=T,
    HF_exp_list=HF,
    HF_exp_err_list=HF_error,
    output_path="/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case/Plots/MAP_value_validation.jpeg",
    std_pred=None,          # or your global std if method="global"
    method="local",
    annotate=True
)





        

        