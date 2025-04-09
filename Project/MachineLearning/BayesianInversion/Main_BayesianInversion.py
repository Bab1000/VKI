import numpy as np
import pandas as pd

from Library.BayesianInversion_vec import *
from Library.Utils_BayesianInversion import *

# ========================================================================
 
# ---------------------------------
# | INPUTS FOR THE BAYESIAN MODEL |
# ---------------------------------

# Bayesian inversion steps
draws = 1100

# Warm-up steps
tune = int(0.20*draws)

# Number of computation chains
chain = 10

# Number of core used
cores = 10

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

# Path previous simulation for restart purposes
# ---------------------------------------------

# If restart is not required --> Leave this spote EMPTY
#restart_sim = f"Results/Posteriors/PostDist_d=1000_c=10_Pstat=4998.226157.nc"

# Path for saving the bayesian inversion simulation
# -------------------------------------------------

# If saving is not required --> Leave this spote EMPTY
save_path = f"Results/Posteriors/PostDist_d={draws}_c={chain}_Pstat={Pstat}.nc"


# Gathering data for bayesian inversion
# -------------------------------------
target_pressure = [5000]
tolerance = 300

Pdyn_test = []
Pstat_test = []
T_test = []
HF_test = []
off_set_Pdyn_test = []
off_set_HF_test = []

for targ_p in target_pressure:
    for i,pstat in enumerate(Pstat):

        if pstat >= targ_p - tolerance and pstat <= targ_p + tolerance:
            Pdyn_test.append(Pdyn[i])
            Pstat_test.append(Pstat[i])
            T_test.append(temperature[i])
            HF_test.append(HF[i])
            off_set_Pdyn_test.append(off_set_Pdyn[i])
            off_set_HF_test.append(off_set_HF[i])


# Initialization of the Bayesian model
model = BayesianInversion(HF_test,sm_q)

for i in range(2):
    # Name pattern
    name_pattern = f"_TC={i}"
    
    model.add_priors("Pdyn" + name_pattern,"normal",Pdyn_test[i],off_set_Pdyn_test[i])
    model.add_priors("Pstat" + name_pattern,"normal",Pstat_test[i],0.05*Pstat_test[i])
    model.add_priors("T"  + name_pattern,"normal",T_test[i],0.05*T_test[i])
    model.add_observed("HF" + name_pattern, "normal", HF_test[i], off_set_HF[i])

model.add_priors("GN","uniform",-4,0)
model.add_priors("GO","uniform",-4,0)

# Running 
trace = model.run_inference(draws,tune,chain,cores)

# Plotting
name_image = f"Results/PostDist_d={draws}_c={chain}_Pstat={Pstat[0]}.jpeg"
#model.plot_posteriors()
model.plot_posteriors_custom(name_image=name_image)

model.predict_from_posterior()

