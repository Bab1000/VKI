import numpy as np
import pandas as pd
import arviz as az  # type: ignore

from Library.BayesianInversion_vec import *
from Library.Utils_BayesianInversion import *


# ========================================================================
 
# ---------------------------------
# | INPUTS FOR THE BAYESIAN MODEL |
# ---------------------------------

# Bayesian inversion steps
draws = 50000

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

# ========================================================================
 
# ----------------------
# | BAYESIAN INVERSION |
# ----------------------

Pdyn = np.linspace(40, 200, 5)               
Pstat = np.full(5, 5000)                      
T = np.linspace(6000, 8000, 5)               
GN = np.full(5, -2)                      
GO = np.full(5, -3)                        

X_vec = np.stack((Pdyn, Pstat, T, GN, GO), axis=1)

HF = sm_q.predict_values(X_vec)[:, 0]  # shape (5,)

# Initialization of the Bayesian model
model = BayesianInversion(HF,sm_q)

for i in range(len(Pdyn)):
    # Name pattern
    name_pattern = f"_TC={i}"
    
    model.add_priors("Pdyn" + name_pattern,"normal",mu = Pdyn[i],uncertainty=0.02*Pdyn[i])
    model.add_priors("Pstat" + name_pattern,"normal",mu=Pstat[i],uncertainty=0.02*Pstat[i])
    model.add_priors("T" + name_pattern,"normal",mu=T[i],uncertainty=0.02*T[i])
    model.add_observed("HF" + name_pattern, "normal", mu=HF[i], uncertainty=0.02*HF[i])

model.add_priors("GN", "uniform", lower=-4, upper=0)
model.add_priors("GO", "uniform", lower=-4, upper=0)


# Global path
# -----------
global_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case_50mbar"

# Path for saving the bayesian inversion simulation
# -------------------------------------------------
# If saving is not required --> Leave this spot EMPTY
save_path_posteriors = os.path.join(global_path,"Posteriors")

# Running 
trace, MAP_values, R_hat_all = model.run_inference(draws,tune,chain,cores,save_path=save_path_posteriors,restart=False)

# Plotting
save_path_results = os.path.join(global_path,"Results")
extension = "jpeg"

model.plot_posteriors_custom(save_path_results,extension)

# MAP/R_hat saving in CSV
SaveAnalysis(save_path_results,MAP_values,R_hat_all)

#model.predict_from_posterior()

# ========================================================================