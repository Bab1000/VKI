import numpy as np
import ctypes
import os
import subprocess
from colorama import Fore, Style, init
from scipy.optimize import minimize_scalar
import pandas as pd
import matplotlib.pyplot as plt

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

        dx_range = np.linspace(1e-6, 1e1, 100)
        qw_vals = []

        """
        for dx in dx_range:
            wdot = (ctypes.c_double * n_species)()
            qw = ctypes.c_double()
            lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx, wdot, ctypes.byref(qw))
            qw_vals.append(qw.value / 1000)

        plt.plot(np.log10(dx_range), np.log(qw_vals))
        plt.axhline(y=np.log(qexp/1000), color='red', linestyle='--', label='qexp')
        plt.xlabel("log10(dx)")
        plt.ylabel("qw [kW/m²]")
        plt.title("qw vs log10(dx)")
        plt.grid(True)
        plt.legend()
        plt.show()
        """

        iteration_count = {"n": 0}
        def log10_equation(log_dx):
            dx = 10 ** log_dx
            if dx <= 0:
                return 1e20  # sécurité : dx doit être positif

            iteration_count["n"] += 1

            # Appel au modèle C++
            wdot = (ctypes.c_double * n_species)()
            qw = ctypes.c_double()
            lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx, wdot, ctypes.byref(qw))

            # Calcul de l’erreur
            qw_kW = qw.value / 1000
            error = abs(qexp / 1000 - qw_kW)

            #print(f"[{iteration_count['n']:02d}] log10(dx) = {log_dx:.3f}, dx = {dx:.3e}, qw = {qw_kW:.3f}, error = {error:.3f}")
            return error

        # Lancement de l'optimisation sur log10(dx)
        result = minimize_scalar(
            log10_equation,
            bounds=(-6, 3),  # dx entre 1e-6 et 1e-1
            method='bounded',
            options={'xatol': 1e-9}
        )

        # Extraction du résultat
        if result.success:
            log_dx_solution = result.x
            solution_dx = 10 ** log_dx_solution
        else:
            print(Fore.RED + "[ERROR] Optimization failed.")
            solution_dx = None

        print(Fore.WHITE + f"---> [INFO] Total iterations = {iteration_count['n']}")

        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()
        lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, solution_dx, wdot, ctypes.byref(qw))

        wdot_np = np.ctypeslib.as_array(wdot, shape=(n_species,))
        print(Fore.GREEN + f"---> [SUCCESS] Results successfully gathered !")

        qw = qw.value

        print(Fore.MAGENTA + f"[RESULTS] Final dx to match experimental qw : {solution_dx}")
        print(Fore.MAGENTA + f"[RESULTS] Final qw (0DReactor) : {qw/1000}")
        print(Fore.MAGENTA + f"[RESULTS] Error : {abs(qexp/1000 - qw/1000)}")

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