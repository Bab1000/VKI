import numpy as np
import pandas as pd
from Utils_SMValidation import *

# ===========================================================================================

# --------------------------
# | LOAD PRE-TRAINED MODEL |
# --------------------------

model_path = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Models/SM_Gamma_Log"

sm_q = LoadModel(model_path)

# ===========================================================================================

# ----------------------------
# | LOAD DATA FROM CSV FILES |
# ----------------------------

StagInput_path = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Validation_Enrico/Data/Data_Enrico_Inputs.xlsx"

StagRes_path = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Validation_Enrico/Data/Data_Enrico_Results.xlsx"


df_FullData = LoadDataEnrico(StagInput_path,StagRes_path)

# ===========================================================================================

# --------------------------------------
# | RUN STAGLINE SM ON VALIDATION DATA |
# --------------------------------------

df_Results = RunSM(sm_q,df_FullData)

# ===========================================================================================

# --------------------
# | PLOTTING RESULTS |
# --------------------

target_pressure = [15, 30, 50, 100, 200]

print(Fore.BLUE + "[STEP] Plotting Results")

for Pc in target_pressure:
    Plot(Pc,df_Results)

PerformanceEval(df_Results)

print(Fore.GREEN + "---> [SUCCESS] Results are successfully plotted in the folder : /home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Data/Validation_Enrico/Results")







