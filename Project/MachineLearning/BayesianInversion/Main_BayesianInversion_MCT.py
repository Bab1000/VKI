import numpy as np
import pandas as pd
import multiprocessing as mlt
from colorama import Fore, Style, init

from Library.BayesianInversion_vec import *
from Library.Utils_BayesianInversion import *

# Initialize colorama
init(autoreset=True)

# ========================================================================
 
# ---------------------------------
# | INPUTS FOR THE BAYESIAN MODEL |
# ---------------------------------

# Bayesian inversion steps
draws = 500000

# Warm-up steps
tune = int(0.20*draws)

# Number of computation chains
chain = 10

# Number of core used
cores = mlt.cpu_count()
print(Fore.YELLOW + f"[CONFIG-INFO] Currently using {cores} cores !")

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
target_pressure = [1500, 5000, 10000, 20000]
target_massflow = [20]
tolerance_press = 300
tolerance_mf = 3

for targ_p in target_pressure:

    for targ_mf in target_massflow:

        # Gathering data for bayesian inversion
        # -------------------------------------
        Pdyn_test = []
        Pstat_test = []
        T_test = []
        HF_test = []
        off_set_Pdyn_test = []
        off_set_HF_test = []

        for i,pstat in enumerate(Pstat):

            if pstat >= targ_p - tolerance_press and pstat <= targ_p + tolerance_press:
                
                if massflow[i] >= targ_mf - tolerance_mf and massflow[i] <= targ_mf + tolerance_mf:

                    Pdyn_test.append(Pdyn[i])
                    Pstat_test.append(Pstat[i])
                    T_test.append(temperature[i])
                    HF_test.append(HF[i])
                    off_set_Pdyn_test.append(off_set_Pdyn[i])
                    off_set_HF_test.append(off_set_HF[i])


        # Initialization of the Bayesian model
        model = BayesianInversion(HF_test,sm_q)

        for i in range(len(Pdyn_test)):
            # Name pattern
            name_pattern = f"_TC={i}"
            
            model.add_priors("Pdyn" + name_pattern,"normal",mu = Pdyn_test[i],uncertainty=off_set_Pdyn_test[i])
            model.add_priors("Pstat" + name_pattern,"normal",mu=Pstat_test[i],uncertainty=0.02*Pstat_test[i])
            model.add_priors("T" + name_pattern,"normal",mu=T_test[i],uncertainty=0.02*T_test[i])
            #model.add_constant("Pdyn" + name_pattern,Pdyn_test[i])
            #model.add_constant("Pstat" + name_pattern,Pstat_test[i])
            #model.add_constant("T" + name_pattern,T_test[i])
            model.add_observed("HF" + name_pattern, "normal", mu=HF_test[i], uncertainty=off_set_HF[i])

        model.add_priors("GN","uniform", lower = -4, upper = 0)
        model.add_priors("GO","uniform", lower = -4, upper = 0)

        # Global path
        # -----------
        global_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Mass_Flow_based_testing/UniformG_MultiCond_Pstat={targ_p/100}/MassFlow={targ_mf}"

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

