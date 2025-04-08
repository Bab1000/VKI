from GA_Functions import *
from Utils_GA_Gamma import *
from colorama import Fore, Style, init
import pdb

import numpy as np

# Initialize colorama for colored terminal output
init(autoreset=True)

# ===========================================================================================

# --------------------------
# | LOAD PRE-TRAINED MODEL |
# --------------------------

model_path = "/home/jpe/VKI/Project/MachineLearning/GeneticAlgorithm/KRGModel/SM_heatFlux"

sm_q = LoadModel(model_path)

# ===========================================================================================

# ---------------------------
# | GA MODEL INITIALIZATION |
# ---------------------------

print(Fore.WHITE + f"---> [INFO] Initialisation of the genetic algorythm model")

# Select a population of n_proc*100 elements
n_p = 100 # Total Population

# Boundaries of the function to be optimized
X_Bounds=[(0,1),(0,1)] 

# Number of iterations
N_ITER=300 

# Proportion of the population initialized using a Gaussian distribution
n_G=0.5

# Interval ratio to calculate standard deviation for Gaussian initialization
sigma_I_r=6

# Initial and Final Mutation Rates 
mu_I=0.02; mu_F=0.0001  

# Portion of the Chromosome subject to Mutation
p_M=0.5 

# Portion subjet to Elitism
n_E = 0.05

# Creating the function to analyze with GA
print(Fore.WHITE + f"---> [INFO] Building function to analyze")
Func = FuncBuild

# ===========================================================================================

# ---------------------------
# | EXPERIMENTAL CONDITIONS |
# ---------------------------

# Static pressure (chamber) [mbar]
target_Pc = [15, 50, 100]

results = []

for target_press in target_Pc:

    print(Fore.WHITE + f"---> [INFO] Looking for data for Pc = {target_press} [Pa]")

    # Path to the CSV
    CSV_path = "/home/jpe/VKI/Project/TestPlasmatron/MassFlowAnalysis/tests_Justin.xlsx"

    exp_HF_data, T_mdot_data, pitot_mdot_data, target_mdot_values = ExpDataLoading(CSV_path,target_press)

    #print(exp_HF_data)
    #pdb.set_trace()

    # ===========================================================================================

    # ----------------
    # | RUN GA MODEL |
    # ----------------

    for mdot_i in target_mdot_values:
        for i in range (len(T_mdot_data[mdot_i])):

            print(Fore.BLUE + "[STEP] Initialization of experimental conditions")

            Pdyn = pitot_mdot_data[mdot_i][i]

            Tinlet = T_mdot_data[mdot_i][i]
            
            Q_target = round(exp_HF_data[mdot_i][i], 3)

            # Run the Animation
            Pc = target_press * 100 # [Pa]

            print(Fore.WHITE + f"---> [INFO] Running model on the following inputs: ")
            print(Fore.WHITE + f"            - mdot = {mdot_i} [g/s]")
            print(Fore.WHITE + f"            - Pc = {Pc} [Pa]")
            print(Fore.WHITE + f"            - Pdyn = {Pdyn} [Pa]")
            print(Fore.WHITE + f"            - Tinlet = {Tinlet} [K]")
            print(Fore.WHITE + f"            - Q_target = {Q_target} [kW]")
         
            X_S, X_U, X_V, Q_pred, Error = GA_model(sm_q,Pdyn,Pc,Tinlet,Func,X_Bounds,Q_target,n_p,N_ITER,n_G,sigma_I_r,mu_I,mu_F,p_M,n_E)

            # Extract gammas
            GammaN = X_S[0]
            GammaO = X_S[1]

            # Store result
            results.append({
                "mdot": mdot_i,
                "Pdyn": Pdyn,
                "Pc": Pc,
                "Tin": Tinlet,
                "gammaN": GammaN,
                "gammaO": GammaO,
                "Q_pred": Q_pred.item(),
                "Q_exp": Q_target,
                "Error": Error
            })
try:
    print(Fore.BLUE + "[STEP] Writting results in CSV")
    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(results)
    path_results = "GA_results_exp_campaign.csv"
    df.to_csv(path_results, index=False, sep=" ")
    print(Fore.GREEN + f"---> [SUCCESS] CSV successfully created : {path_results}")
except Exception as e:
    print(Fore.GREEN + f"---> [ERROR] Error while writting results : {e}")



# ===========================================================================================

# -------------------------------
# | RUN GA MODEL WITH ANIMATION |
# -------------------------------

#Q_target = 1547

# Run the Animation
#X_S, X_U, X_V=Anim_COMP(sm_q,Pdyn,Pc,Tinlet,Func,X_Bounds,Q_target,n_p,N_ITER,n_G,sigma_I_r,mu_I,mu_F,p_M,n_E,x_1m=0,x_1M=1,x_2m=0,x_2M=1,npoints=100,Name_Video='GA_gamma.gif')