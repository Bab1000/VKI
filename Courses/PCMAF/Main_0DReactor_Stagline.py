from Utils_0DReactor import *
import matplotlib.pyplot as plt
import pandas as pd
import os


# === Compilation + Loading of the model ===

cpp_file = "/home/jpe/VKI/Courses/PCMAF/0DReactor/0DReactor.cpp"
exe_file = "0DReactor.so"
mutationpp_include = "/home/jpe/Mutationpp/src"
eigen_include = "/home/jpe/Mutationpp/thirdparty/eigen"
lib_path = "/home/jpe/Mutationpp/build/src"
lib_name = "-lmutation++"

Load0DReactor(cpp_file, exe_file, mutationpp_include, eigen_include, lib_path, lib_name)

# ===========================================================================================

# ----------------
# | DATA LOADING |
# ----------------

# File paths for training and validation data
data_folder_name = "/home/jpe/VKI/Courses/PCMAF/Data"
train_Y_file = data_folder_name + "/YTraining.dat"
train_X_file = data_folder_name + "/Training_Xdata.dat"
valid_Y_file = data_folder_name + "/YValidation.dat"
valid_X_file = data_folder_name + "/Validation_Xdata.dat"


# Data containers
XT, YT, XV, YV = [], [], [], []

# Load training data
load_data(train_Y_file, train_X_file, XT, YT, desc="Loading Training Data")

print(Fore.GREEN + "[SUCCESS] Training data successfully loaded!")

# ===========================================================================================

# === Runing 0DReactor ===

results_path = "/home/jpe/VKI/Courses/PCMAF/Results"

path_csv = os.path.join(results_path,f"Res_tuning_Pstat_Stagline_Single_dx.csv")

Results = []
inputs = []
HF = []

for idx, line in enumerate(XT): 
    if line[2] > 7000:
        inputs.append(line)
        HF.append(YT[idx])

for i in range(len(inputs)):
    input = inputs[i]
    mixture = f"air_11"
    Tinlet = input[2]       # [K]
    Tsurface = 350.0        # [K]
    Pstat_in = input[1]     # [Pa]
    n_species = 11          # [-]
    qexp = HF[i]            # [w/m2]
    gamma_n = 10**input[3]
    gamma_o = 10**input[4]

    update_gsi_gamma("gsi_surface_cat.xml",gamma_n,gamma_o)

    wdot, qw, dx_diff, dx_conv, error = run0DReactorSingledx(mixture,Tinlet,Tsurface,Pstat_in,qexp,exe_file,n_species)

    if error < qexp /1000 * 0.2:
        Results.append({
            "index": i,
            "Pstat": Pstat_in,
            "Tinlet": Tinlet,
            "Tsurface": Tsurface,
            "qexp": qexp,
            "dx_diff": dx_diff,
            "dx_conv": dx_conv,
            "qw": qw
        })
    else:
        print(Fore.YELLOW + "[RESULTS-INFO] No dx found giving errors below 20% ")


df = pd.DataFrame(Results)

df.to_csv(path_csv, index=False)



        