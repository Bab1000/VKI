import numpy as np
import ctypes
import os
import subprocess
from colorama import Fore, Style, init
from scipy.optimize import newton
import pandas as pd

# Initialize colorama
init(autoreset=True)

def Load0DReactor(cpp_file, exe_file, mutationpp_include, eigen_include, lib_path, lib_name):
    # === COMPILATION ===
    print(Fore.BLUE + f"[STEP] Compilation of the model : {exe_file}")

    print(Fore.WHITE + f"---> [INFO] Preparing the libraries ...")

    compile_command = [
        "g++", "-shared", "-fPIC", cpp_file,
        "-o", exe_file,
        "-I" + os.path.join(mutationpp_include, "thermo"),
        "-I" + os.path.join(mutationpp_include, "general"),
        "-I" + eigen_include,
        "-I" + os.path.join(mutationpp_include, "kinetics"),
        "-I" + os.path.join(mutationpp_include, "utilities"),
        "-I" + os.path.join(mutationpp_include, "numerics"),
        "-I" + os.path.join(mutationpp_include, "transport"),
        "-I" + os.path.join(mutationpp_include, "transfer"),
        "-I" + os.path.join(mutationpp_include, "gsi"),
        "-L" + lib_path, lib_name
    ]

    print(Fore.WHITE + f"---> [INFO] Compilation ...")
    compilation = subprocess.run(compile_command, capture_output=True, text=True)

    if compilation.returncode != 0:
        print(Fore.RED + "[ERROR] Error during the compilation !")
        print(Fore.RED + "STDOUT:\n", compilation.stdout)
        print(Fore.RED + "STDERR:\n", compilation.stderr)
        exit(1)

    print(Fore.GREEN + f"---> [SUCCESS] Compilation successfully executed !")

def run0DReactor(mixture,Tinlet,Tsurface,Pstat,qexp,exe_file,n_species):

    try:

        print(Fore.BLUE + f"[STEP] Running the model : {exe_file}")

        print(Fore.WHITE + f"---> [INFO] Loading the model ...")
        lib = ctypes.CDLL("./" + exe_file)
        lib.ZeroDReactor.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
        lib.ZeroDReactor.restype = ctypes.c_int
        print(Fore.GREEN + f"---> [SUCCESS] Model successfully loaded !")

        mixture = mixture.encode("utf-8")

        print(Fore.WHITE + f"---> [INFO] Defining Newton function ...")

        iteration_count = {"n": 0}
        def equation(dx):
            iteration_count["n"] += 1
            wdot = (ctypes.c_double * n_species)()
            qw = ctypes.c_double()
            lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx, wdot, ctypes.byref(qw))
            return qexp - qw.value

        dx_initial_guess = 0.001
        solution_dx = newton(equation, dx_initial_guess, tol=1e-8)

        print(Fore.WHITE + f"---> [INFO] Final call with solution dx = {solution_dx:.6e}")
        print(Fore.WHITE + f"---> [INFO] Newton iterations = {iteration_count['n']}")

        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()
        lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, solution_dx, wdot, ctypes.byref(qw))

        wdot_np = np.ctypeslib.as_array(wdot, shape=(n_species,))
        print(Fore.GREEN + f"---> [SUCCESS] Results successfully gathered !")

        qw = qw.value

        print(Fore.MAGENTA + f"[RESULTS] Final dx to match experimental qw : {solution_dx}")
        print(Fore.MAGENTA + f"[RESULTS] Final qw (0DReactor) : {qw/1000}")

        return wdot_np, qw, solution_dx

    except Exception as e:
        print(Fore.RED + f"[ERROR] Error while running {exe_file}: {e}")
        return None, None, None
    
def CSVReader(CSV_path,columns_to_check):
    try:

        print(Fore.BLUE + "[STEP] Loading data from CSV")

        # Reading the CSV file with all the data
        df = pd.read_excel(CSV_path, engine="openpyxl")

        print(Fore.WHITE + "---> [INFO] Treating data ...")
        # Replace "NA" strings with actual NaN values
        df.replace("NA", np.nan, inplace=True)

        # Drop rows where any of the specified columns contain NaN
        df_cleaned = df.dropna(subset=columns_to_check)

        print(Fore.WHITE + "---> [INFO] Preparing the data for external use ...")
        # Gathering the cleaned data
        pressure = np.array(df_cleaned["Pressure[mbar]"].dropna()) * 100 # [mbar] --> [Pa]
        massflow = np.array(df_cleaned["massflow [g/s]"].dropna())
        power = np.array(df_cleaned["Power[kW]"].dropna())
        heat_flux = np.array(df_cleaned["HeatFlux(HS50mm)[kW/m2]"].dropna()) * 1000 # [kW/m2] --> [W/m2]
        off_set_heat_flux = np.array(df_cleaned["offsetHF[kW/m2]"].dropna()) * 1000 # [kW/m2] --> [W/m2]
        pitot = np.array(df_cleaned["Pitot[Pa]"].dropna())
        off_set_pitot = np.array(df_cleaned["offsetPitot[Pa]"].dropna())
        temperature = np.array(df_cleaned["T [K] (x = 375mm, r = 0mm)"].dropna())

        # Ensure + values in off_sets
        off_set_heat_flux = abs(off_set_heat_flux)
        off_set_pitot = abs(off_set_pitot)

        print(Fore.GREEN + "---> [SUCCESS] Data loaded successfully !")

        return pressure,massflow,power,heat_flux,off_set_heat_flux,pitot,off_set_pitot,temperature
    
    except Exception as e:
        print(Fore.RED + f"[ERROR] Data loading failed: {e}")
