import numpy as np
import ctypes
import os
import subprocess
from colorama import Fore, Style, init
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import pandas as pd
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm
import xml.etree.ElementTree as ET

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

def run0DReactorMultidx(mixture,Tinlet,Tsurface,Pstat,qexp,exe_file,n_species):

    print(Fore.BLUE + f"[STEP] Running the model : {exe_file}")

    print(Fore.WHITE + f"---> [INFO] Loading the model ...")
    lib = ctypes.CDLL("./" + exe_file)
    lib.ZeroDReactor.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.ZeroDReactor.restype = ctypes.c_int
    print(Fore.GREEN + f"---> [SUCCESS] Model successfully loaded !")

    mixture = mixture.encode("utf-8")

    print(Fore.WHITE + f"---> [INFO] Defining Newton function ...")

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

    # Initialisation
    iteration_count = {"n": 0, "best": None}

    # Error equation
    def log10_equation(log_dx_diff, log_dx_conv):

        # Duplacating dx for convection and diffusion
        dx_diff = 10 ** log_dx_diff
        dx_conv = 10 ** log_dx_conv

        # Penalising negative values for both dx
        if dx_diff <= 0 or dx_conv <= 0:
            return 1e20

        iteration_count["n"] += 1
        
        # Preparing storage for wdot and qw
        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()

        # Calling 0Drecator in C++
        try:
            lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx_diff, dx_conv, wdot, ctypes.byref(qw))
        except Exception as e:
            print(Fore.RED + f"[ERROR] Call to ZeroDReactor failed: {e}")
            return 1e20

        # Working with kW
        qw_kW = qw.value / 1000
        error = abs(qexp / 1000 - qw_kW)

        #print(Fore.YELLOW + f"---> [TUNING] Iter {iteration_count['n']:03d} | dx_diff = {dx_diff:.3e} | dx_conv = {dx_conv:.3e} | error = {error:.5f}")

        # Keeping a trace of the best value for dx_conv and dx_diff
        if iteration_count["best"] is None or error < iteration_count["best"]["error"]:
            iteration_count["best"] = {
                "dx_diff": dx_diff,
                "dx_conv": dx_conv,
                "error": error
            }

        # Computing error
        if error < qexp /1000 *0.2:  # Stops if error is lower than 20%
            raise EarlyStop()

        return error

    # Wrapper for the main error function
    def log10_equation_wrapper(log_dx):
        return log10_equation(*log_dx)

    # Initialisation of the minimisation process
    initial_guess = [-3, -3]
    bounds = [(-4, 0), (-4, 0)]

    # Minimization
    try:
        result = minimize(
            log10_equation_wrapper,
            x0=initial_guess,
            bounds=bounds,
            method='Powell',
            options={'ftol': 1e-5}
        )
    except EarlyStop:
        print(Fore.CYAN + "[INFO] Early stopping: error below threshold reached.")
        result = None


    # Get best point
    if iteration_count["best"]:
        solution_dx_diff = iteration_count["best"]["dx_diff"]
        solution_dx_conv = iteration_count["best"]["dx_conv"]
        final_error = iteration_count["best"]["error"]

        print(Fore.WHITE + f"---> [INFO] Total iterations = {iteration_count['n']}")
        print(Fore.GREEN + f"---> [SUCCESS] Best error = {final_error:.5f}")

        # Final evaluation of qw predicted with new dx
        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()
        lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, solution_dx_diff, solution_dx_conv, wdot, ctypes.byref(qw))
        wdot_np = np.ctypeslib.as_array(wdot, shape=(n_species,))
        qw = qw.value

        print(Fore.MAGENTA + f"[RESULTS] Final dx_diff to match experimental qw : {solution_dx_diff}")
        print(Fore.MAGENTA + f"[RESULTS] Final dx_conv to match experimental qw : {solution_dx_conv}")
        print(Fore.MAGENTA + f"[RESULTS] Final qw (0DReactor) : {qw / 1000}")
        print(Fore.MAGENTA + f"[RESULTS] Error : {abs(qexp / 1000 - qw / 1000)}")

        # Working in kW
        qw_kW = qw / 1000
        error = abs(qexp / 1000 - qw_kW)

        return wdot_np, qw, solution_dx_diff, solution_dx_conv,error

    else:
        print(Fore.RED + "[ERROR] No valid dx values found.")
        return None, None, None, None

def run0DReactorSingledx(mixture,Tinlet,Tsurface,Pstat,qexp,exe_file,n_species):

    print(Fore.BLUE + f"[STEP] Running the model : {exe_file}")

    print(Fore.WHITE + f"---> [INFO] Loading the model ...")
    lib = ctypes.CDLL("./" + exe_file)
    lib.ZeroDReactor.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.ZeroDReactor.restype = ctypes.c_int
    print(Fore.GREEN + f"---> [SUCCESS] Model successfully loaded !")

    mixture = mixture.encode("utf-8")

    print(Fore.WHITE + f"---> [INFO] Defining Newton function ...")

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

    # Initialisation
    iteration_count = {"n": 0, "best": None}

    # main error function
    def log10_equation(log_dx):  
        dx_diff = 10 ** log_dx
        dx_conv = 10 ** log_dx

        # Penalizing negative values for dx
        if dx_diff <= 0 or dx_conv <= 0:
            return 1e20

        iteration_count["n"] += 1

        # Preparing storage
        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()

        # Ruining 0D reactor in C++
        try:
            lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx_diff, dx_conv, wdot, ctypes.byref(qw))
        except Exception as e:
            print(Fore.RED + f"[ERROR] Call to ZeroDReactor failed: {e}")
            return 1e20

        # Working in kW
        qw_kW = qw.value / 1000
        error = abs(qexp / 1000 - qw_kW)

        #print(Fore.YELLOW + f"---> [TUNING] Iter {iteration_count['n']:03d} | dx_diff = {dx_diff:.3e} | dx_conv = {dx_conv:.3e} | error = {error:.5f}")

        # Keeping trace of the best dx
        if iteration_count["best"] is None or error < iteration_count["best"]["error"]:
            iteration_count["best"] = {
                "dx_diff": dx_diff,
                "dx_conv": dx_conv,
                "error": error
            }   

        # Computing errors
        if error < qexp /1000 * 0.2:  # Stops when error is lower than 20%
            raise EarlyStop()

        return error

    # Error minimization
    try:
        result = minimize_scalar(
            log10_equation,
            bounds=(-6, 0),
            method='bounded',
            options={'xatol': 1e-5}
        )
    except EarlyStop:
        print(Fore.CYAN + "[INFO] Early stopping: error below threshold reached.")
        result = None


    # Get best dx
    if iteration_count["best"]:
        solution_dx_diff = iteration_count["best"]["dx_diff"]
        solution_dx_conv = iteration_count["best"]["dx_conv"]
        final_error = iteration_count["best"]["error"]

        print(Fore.WHITE + f"---> [INFO] Total iterations = {iteration_count['n']}")
        print(Fore.GREEN + f"---> [SUCCESS] Best error = {final_error:.5f}")

        # Final evaluation of qw predicted with new dx
        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()
        lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, solution_dx_diff, solution_dx_conv, wdot, ctypes.byref(qw))
        wdot_np = np.ctypeslib.as_array(wdot, shape=(n_species,))
        qw = qw.value

        print(Fore.MAGENTA + f"[RESULTS] Final dx_diff to match experimental qw : {solution_dx_diff}")
        print(Fore.MAGENTA + f"[RESULTS] Final dx_conv to match experimental qw : {solution_dx_conv}")
        print(Fore.MAGENTA + f"[RESULTS] Final qw (0DReactor) : {qw / 1000}")
        print(Fore.MAGENTA + f"[RESULTS] Error : {abs(qexp / 1000 - qw / 1000)}")

        # Working in kW
        qw_kW = qw / 1000
        error = abs(qexp / 1000 - qw_kW)

        return wdot_np, qw, solution_dx_diff, solution_dx_conv, error

    else:
        print(Fore.RED + "[ERROR] No valid dx values found.")
        return None, None, None, None, None

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

def load_data(y_file, x_file, X_list, Y_list, desc="Loading Data"):
    """Loads input and output data from files while showing a tqdm progress bar."""
    try:
        with open(y_file) as y_f, open(x_file) as x_f:
            y_lines = y_f.readlines()
            x_lines = x_f.readlines()

            for line_index, line in tqdm(enumerate(y_lines), total=len(y_lines), desc=desc):
                if line_index > 0:
                    splitted = line.split()
                    
                    # Handling non-converged cases
                    if splitted[2] == "NotConv":
                        splitted[2] = 999
                    
                    conv = float(splitted[2])
                    
                    # If convergence is acceptable, add the data
                    if conv < 1e-2:
                        Y_list.append(float(splitted[1]))
                        x_values = x_lines[line_index].split()
                        X_list.append([float(x_values[1]), float(x_values[2]), float(x_values[3]), np.log10(float(x_values[4])), np.log10(float(x_values[5]))])

        print(Fore.GREEN + "---> [SUCCESS] Data loaded successfully!")

    except FileNotFoundError as e:
        print(Fore.RED + f"---> [ERROR] File not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Failed to load validation data: {e}")
        sys.exit(1)

def update_gsi_gamma(xml_file_path, gamma_O, gamma_N, output_file_path=None):
    """
    Updates the gamma_const values for O and N in a GSI XML file.
    """

    # Reading xml file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    updated = {"O": False, "N": False}

    # Update gamma values
    for reaction in root.findall(".//reaction"):
        gamma_elem = reaction.find("gamma_const")
        if gamma_elem is not None:
            gamma_text = gamma_elem.text.strip()
            if "O:" in gamma_text:
                gamma_elem.text = f"O:{gamma_O:.5f}"
                updated["O"] = True
            elif "N:" in gamma_text:
                gamma_elem.text = f"N:{gamma_N:.5f}"
                updated["N"] = True

    # Save XML first
    if output_file_path is None:
        output_file_path = xml_file_path  # overwrite original

    tree.write(output_file_path, encoding="utf-8", xml_declaration=False)

    # Then fix &gt; to => 
    with open(output_file_path, "r", encoding="utf-8") as f:
        content = f.read()

    content = content.replace("=&gt;", "=>").replace("&gt;", ">")

    with open(output_file_path, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"---> [INFO] XML updated and saved to: {output_file_path}")
    print(f"---> [DEBUG] Updated gamma_O: {updated['O']}, gamma_N: {updated['N']}")

def run0DReactorSingle_dxConvOnly(mixture, Tinlet, Tsurface, Pstat, qexp, exe_file, n_species):

    print(Fore.BLUE + f"[STEP] Running the model : {exe_file}")

    print(Fore.WHITE + f"---> [INFO] Loading the model ...")
    lib = ctypes.CDLL("./" + exe_file)
    lib.ZeroDReactor.argtypes = [
        ctypes.c_char_p, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
    ]
    lib.ZeroDReactor.restype = ctypes.c_int
    print(Fore.GREEN + f"---> [SUCCESS] Model successfully loaded !")

    mixture = mixture.encode("utf-8")

    print(Fore.WHITE + f"---> [INFO] Defining Newton function ...")

    # === Initialization ===
    iteration_count = {"n": 0, "best": None}
    dx_diff_fixed = 1e-3

    # === Error function with fixed dx_diff ===
    def log10_equation(log_dx_conv):
        dx_conv = 10 ** log_dx_conv

        if dx_conv <= 0:
            return 1e20

        iteration_count["n"] += 1

        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()

        try:
            lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx_diff_fixed, dx_conv, wdot, ctypes.byref(qw))
        except Exception as e:
            print(Fore.RED + f"[ERROR] Call to ZeroDReactor failed: {e}")
            return 1e20

        qw_kW = qw.value / 1000
        error = abs(qexp / 1000 - qw_kW)

        # Update best solution if better
        if iteration_count["best"] is None or error < iteration_count["best"]["error"]:
            iteration_count["best"] = {
                "dx_diff": dx_diff_fixed,
                "dx_conv": dx_conv,
                "error": error
            }

        # Early stop if acceptable
        #if error < qexp / 1000 * 0.2:
        #    raise EarlyStop()

        return error

    # === Optimization ===
    try:
        result = minimize_scalar(
            log10_equation,
            bounds=(-6, 0),
            method='bounded',
            options={'xatol': 1e-5}
        )
    except EarlyStop:
        print(Fore.CYAN + "[INFO] Early stopping: error below threshold reached.")
        result = None

    # === Final evaluation ===
    if iteration_count["best"]:
        dx_diff = iteration_count["best"]["dx_diff"]
        dx_conv = iteration_count["best"]["dx_conv"]
        final_error = iteration_count["best"]["error"]

        print(Fore.WHITE + f"---> [INFO] Total iterations = {iteration_count['n']}")
        print(Fore.GREEN + f"---> [SUCCESS] Best error = {final_error:.5f}")

        # Final call to get wdot and qw
        wdot = (ctypes.c_double * n_species)()
        qw = ctypes.c_double()
        lib.ZeroDReactor(mixture, Tinlet, Tsurface, Pstat, dx_diff, dx_conv, wdot, ctypes.byref(qw))
        wdot_np = np.ctypeslib.as_array(wdot, shape=(n_species,))
        qw = qw.value

        print(Fore.MAGENTA + f"[RESULTS] Final dx_diff (fixed) : {dx_diff}")
        print(Fore.MAGENTA + f"[RESULTS] Final dx_conv (optimized) : {dx_conv}")
        print(Fore.MAGENTA + f"[RESULTS] Final qw (0DReactor) : {qw / 1000}")
        print(Fore.MAGENTA + f"[RESULTS] Error : {abs(qexp / 1000 - qw / 1000)}")

        return wdot_np, qw, dx_diff, dx_conv, final_error

    else:
        print(Fore.RED + "[ERROR] No valid dx_conv value found.")
        return None, None, None, None, None


class EarlyStop(Exception):
    pass