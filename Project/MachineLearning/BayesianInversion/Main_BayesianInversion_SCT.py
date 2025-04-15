import numpy as np
import pandas as pd

from Library.BayesianInversion_vec import *
from Library.Utils_BayesianInversion import *

# ========================================================================
 
# ---------------------------------
# | INPUTS FOR THE BAYESIAN MODEL |
# ---------------------------------

# Bayesian inversion steps
draws = 1000

# Warm-up steps
tune = int(0.20*draws)

# Number of computation chains
chain = 5

# Number of core used
cores = 5

# ========================================================================
 
# ------------------------------------
# | STAGLINE SURROGATE MODEL LOADING |
# ------------------------------------

model_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/KRGModel/SM_Gamma_Log"

sm_q = LoadModel(model_path)

# Creating an instance of the surrogate model


# ========================================================================
 
# ----------------
# | DATA LOADING |
# ----------------

# Path to test campaign data
CSV_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Data/tests-mass-flow-Justin.xlsx"

# Columns to check for wrong data
columns_to_check = ["Pressure[mbar]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)", "HeatFlux(HS50mm)[kW/m2]"]

Pstat,massflow,power,HF,off_set_HF,Pdyn,off_set_Pdyn,temperature = CSVReader(CSV_path,columns_to_check)

# ========================================================================
 
# ----------------------
# | BAYESIAN INVERSION |
# ----------------------

# Target static pressure 
# ----------------------
target_pressure = [5000]
target_massflow = 16
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

for targ_p in target_pressure:
    for i,pstat in enumerate(Pstat):

        if pstat >= targ_p - tolerance_press and pstat <= targ_p + tolerance_press:
            
            if massflow[i] >= target_massflow - tolerance_mf and massflow[i] <= target_massflow + tolerance_mf:

                Pdyn_test.append(Pdyn[i])
                Pstat_test.append(Pstat[i])
                T_test.append(temperature[i])
                HF_test.append(HF[i])
                off_set_Pdyn_test.append(off_set_Pdyn[i])
                off_set_HF_test.append(off_set_HF[i])


for i in range(len(Pdyn_test)):
    # Name pattern
    #name_pattern = f"_TC={i}"

    # Initialization of the Bayesian model
    model = BayesianInversion(HF_test[i],sm_q)
    
    model.add_priors("Pdyn","normal",mu = Pdyn_test[i],uncertainty=off_set_Pdyn_test[i])
    model.add_priors("Pstat","normal",mu=Pstat_test[i],uncertainty=0.10*Pstat_test[i])
    model.add_priors("T","normal",mu=T_test[i],uncertainty=0.10*T_test[i])
    model.add_observed("HF", "normal", mu=HF_test[i], uncertainty=off_set_HF[i])

    model.add_priors("GN","uniform", lower = -4, upper = 0)
    model.add_priors("GO","uniform", lower = -4, upper = 0)

    # Global path
    # -----------

    global_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Single_Condition_Testing/"

    # Path previous simulation for restart purposes
    # ---------------------------------------------
    # If restart is not required --> Leave this spot EMPTY
    #restart_sim = global_path + f"Posteriors/PostDist_Unit_d={draws}_c={chain}_Pstat={target_pressure[0]/100}_T={T_test[i]}.nc"

    # Path for saving the bayesian inversion simulation
    # -------------------------------------------------
    # If saving is not required --> Leave this spot EMPTY
    save_path_posteriors = global_path + f"Posteriors/PostDist_Unit_d={draws}_c={chain}_Pstat={target_pressure[0]/100}_T={T_test[i]}.nc"

    # Running 
    trace, MAP_values, R_hat_all = model.run_inference(draws,tune,chain,cores,save_path=save_path_posteriors)

    # Plotting
    save_path_results = global_path + f"Results/PostDist_Unit_d={draws}_c={chain}_Pstat={target_pressure[0]/100}_T={T_test[i]}"
    extension = "jpeg"

    model.plot_posteriors_custom(save_path_results,extension)

    # MAP/R_hat saving in CSV
    SaveAnalysis(save_path_results,MAP_values,R_hat_all)

    #model.predict_from_posterior()

