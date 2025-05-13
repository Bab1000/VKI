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

# === Compilation + Loading of the model ===

# Path to test campaign data
CSV_path = "/home/jpe/VKI/Courses/PCMAF/Data/tests-mass-flow-Justin.xlsx"

# Columns to check for wrong data
columns_to_check = ["Pressure[mbar]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)", "HeatFlux(HS50mm)[kW/m2]"]

Pstat,massflow,power,heat_flux,off_set_heat_flux,pitot,off_set_pitot,temperature = CSVReader(CSV_path,columns_to_check)

# === Runing 0DReactor ===

target_pressure = [1500,5000,10000,20000]
tolerance = 300

print(np.shape(Pstat))

for targ_p in target_pressure:

    Pstat_test = []
    T_test = []
    HF_test = []

    results_path = "/home/jpe/VKI/Courses/PCMAF/Results"

    path_csv = os.path.join(results_path,f"Res_tuning_Pstat={targ_p}.csv")

    Results = []

    for i,pstat in enumerate(Pstat):

        if pstat >= targ_p - tolerance and pstat <= targ_p + tolerance:
            Pstat_test.append(Pstat[i])
            T_test.append(temperature[i])
            HF_test.append(heat_flux[i])

    for i in range(len(Pstat_test)):
        mixture = f"air_11_{int(targ_p/100)}mbar"
        Tinlet = T_test[i]        # [K]
        Tsurface = 350.0          # [K]
        Pstat_in = Pstat_test[i]     # [Pa]
        n_species = 11            # [-]
        qexp = HF_test[i]         # [w/m2]

        wdot, qw, dx = run0DReactor(mixture,Tinlet,Tsurface,Pstat_in,qexp,exe_file,n_species)

        Results.append({
        "Pstat": targ_p,
        "Tinlet": Tinlet,
        "Tsurface": Tsurface,
        "qexp": qexp,
        "dx": dx,
        "qw": qw
    })

    df = pd.DataFrame(Results)

    df.to_csv(path_csv, index=False)



        