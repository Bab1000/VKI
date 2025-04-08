import shutil
import os
import subprocess
import numpy as np
from scipy.optimize import newton
from colorama import Fore, Style, init
import time
import matplotlib.pyplot as plt
import pdb
import glob
import pandas as pd

init(autoreset=True)

def UserInputsProcessing(Txin, pc, mixture, pdyn, uin, vin, R, n_species, cfl_val, Twall, cfl_inter, cfl_adaptive, Log_CFL,residual,restart):
    # Initialisation of the inputs and init_cond dictionary
    # -----------------------------------------------------
    print(Fore.BLUE + "[STEP] Processing user inputs")
    
    inputs = {}
    init_cond = {}

    # Locating the micture file from M++
    # ----------------------------------
    mixture_file_name = mixture + f"_{int(pc/100)}mbar"
    
    # Check for type of inputs (Pdyn or velocity)
    # -------------------------------------------
    if pdyn is not None and uin is None and vin is None:
        print(Fore.WHITE + "   [INFO] Calculating velocity from dynamic pressure")
        
        # Computing the mixture density
        rho = MPPMixtureDensity(Txin, pc, mixture_file_name)
        print(Fore.WHITE + f"   [INFO] Computed density: {rho}")

        # Computing the mixture dynamic viscosity
        mu = MPPViscosity(Txin, pc, mixture_file_name)
        print(Fore.WHITE + f"   [INFO] Computed viscosity: {mu}")

        # Computing the velocity
        U = PdynToVelocity(pdyn, rho, mu, R)
        uin = -U
        vin = U
        init_cond["Velocities"] = np.array([uin, vin])
        print(Fore.WHITE + f"   [INFO] Computed velocity: {U}")
        
    
    elif pdyn is None and uin is not None and vin is not None:
        print(Fore.WHITE + "   [INFO] Using provided velocity values")
        init_cond["Velocities"] = np.array([uin, vin])

    elif pdyn is None and uin is None and vin is None:
        print(Fore.RED + "[ERROR] The value of the inlet dynamic pressure or velocity must be specified")

    elif pdyn is not None and uin is not None and vin is not None:
        print(Fore.RED + "[ERROR] Either the inlet dynamic pressure or velocity must be specified, but not both.")

    else:
        print(Fore.RED + "[ERROR] Wrong value specification for the dynamic pressure or the velocity")

    # Name of the simulation
    # ----------------------
    sim_name = f"sim_Pc={pc}_T={Txin}_U={np.abs(uin):.1f}_mix={mixture_file_name}"
    
    # Processing the rest of the inputs
    # ---------------------------------
    inputs.update({
        "Simulation_Name": sim_name,
        "Mixture": mixture_file_name,
        "Number_of_species": n_species,
        "CFL_number": cfl_val,
        "Twall": Twall,
        "Inter_CFL": cfl_inter,
        "Adaptive_CFL": cfl_adaptive,
        "Log_CFL": Log_CFL,
        "Residual": residual,
        "Restart": restart
    })
    
    # Gathering the species density
    density = MPPSpeciesDensities(Txin, pc, mixture_file_name)
    print(Fore.WHITE + f"   [INFO] Computed species densities: {density}")

    # Gathering the mass fractions
    mass_fractions = MPPMassFractions(Txin, pc, mixture_file_name)
    print(Fore.WHITE + f"   [INFO] Computed mass fractions: {mass_fractions}")
    
    # Initial conditions values
    init_cond["Pressure"] = pc
    init_cond["Densities"] = density
    init_cond["Mass fractions"] = mass_fractions
    init_cond["Temperatures"] = np.array([Txin])
    
    return inputs, init_cond, sim_name, mixture

def InputFileGenerator(stagline_simulations_path,stagline_restart_simulations_path, input_template_path, catalicity_files_path, sim_name, inputs, init_cond,mixture,air_5_restart):
    """
    Generate and modify the input file for a Stagline simulation.
    """

    print(Fore.BLUE + "[STEP] Generation of the input file")

    # Creating simulation directory
    if not os.path.isdir(stagline_simulations_path):
        os.makedirs(stagline_simulations_path)
        print(Fore.WHITE + f"   [INFO] Created simulation directory: {stagline_simulations_path}")
    else:
        print(Fore.WHITE + f"   [INFO] Simulation directory found: {stagline_simulations_path}")
    

    # Name of the simulation folder
    sim_folder_path = os.path.join(stagline_simulations_path, sim_name)

    # Creating simulation folder
    if not os.path.isdir(sim_folder_path):
        os.makedirs(sim_folder_path)
        print(Fore.WHITE + f"   [INFO] Created simulation folder: {sim_folder_path}")
    else:
        print(Fore.YELLOW + f"   [WARNING] Simulation folder already exists: {sim_folder_path}")
    
    # Copying necessary files from macro folder
    
    try:
        shutil.copy(os.path.join(input_template_path, f"Example_input_sub_{mixture}"), sim_folder_path)
        shutil.copy(os.path.join(input_template_path, "mesh.dat"), sim_folder_path)
        shutil.copy(os.path.join(input_template_path, "cfl"), sim_folder_path)
        shutil.copy(os.path.join(input_template_path, "cfl_history_log"), sim_folder_path)
        print(Fore.WHITE + "   [INFO] Template files successfully copied.")
    except FileNotFoundError:
        print(Fore.RED + "[ERROR] One or more template files are missing!")
        return None
    
    # Management of the catalicity files

    # Template for mixture file
    mixture_temp = f"{mixture}_{int(init_cond['Pressure']/100)}mbar.xml"

    # Template for the mechanism_file
    mechanism_temp = f"{mixture}_mech.xml"

    # Template for catalicity file
    catalicity_temp = f"gsi_surface_cat_{int(init_cond['Pressure']/100)}mbar.xml"

    
    try:
        shutil.copy(os.path.join(catalicity_files_path, mixture_temp), sim_folder_path)
        shutil.copy(os.path.join(catalicity_files_path, mechanism_temp), sim_folder_path)
        shutil.copy(os.path.join(catalicity_files_path, catalicity_temp), sim_folder_path)
        print(Fore.WHITE + "   [INFO] Mixture, mechanism and catalicity files successfully copied.")
    except FileNotFoundError:
        print(Fore.RED + "[ERROR] One or more mixture/mechanism/catalicity files are missing!")
        return None

    
    # Renaming the input file
    sim_input = os.path.join(sim_folder_path, f"Example_input_sub_{mixture}")
    renamed_input = os.path.join(sim_folder_path, "input")
    os.rename(sim_input, renamed_input)
    
    print(Fore.WHITE + "   [INFO] Modification of the input file with user inputs")

    with open(renamed_input, "r") as input_file:
        lines = input_file.readlines()
    
    # replacing inputs conditions
    for key, value in inputs.items():
        replace_inputs(lines, key, str(value))
    
    # Replacing initial conditions
    replace_init_cond(lines, init_cond)
    
    with open(renamed_input, 'w') as input_file:
        input_file.writelines(lines)
    
    print(Fore.GREEN + "   [SUCCESS] Input file successfully modified.")

    # Use previous Air_5 converged simulation
    if air_5_restart == ".TRUE.":
        RestartFromPreviousSim(stagline_simulations_path,stagline_restart_simulations_path,sim_name,mixture)


    return sim_folder_path

def RestartFromPreviousSim(stagline_simulations_global_path,stagline_restart_simulations_path, sim_name, mixture):
    # Getting the initial file path of the air_5 simulation
    #old_sim_name = sim_name.replace(mixture, "air_5") 
    #old_sim_folder = os.path.join(stagline_simulations_global_path, old_sim_name)
    #old_sim_folder = old_sim_folder.replace(mixture, "air_5") 

    print(Fore.WHITE + f"   [INFO] Looking for air_5 restart files in: {stagline_restart_simulations_path}")

    # Define the search patterns
    pattern1 = "*restart.dat"
    pattern2 = "*restartBC.dat"

    # Search for files matching the pattern
    file1_list = glob.glob(os.path.join(stagline_restart_simulations_path, pattern1))
    file2_list = glob.glob(os.path.join(stagline_restart_simulations_path, pattern2))

    # Check if any files were found
    if file1_list and file2_list:
        
        # Get the first matching file (assuming only one should be used)
        file1 = file1_list[0]
        file2 = file2_list[0]

        # Copy the restart files into the current simulation folder
        current_sim_folder = os.path.join(stagline_simulations_global_path, sim_name)

        try:
            # Define new filenames by replacing "air_5" with "air_7"
            new_file1_name = sim_name + "_restart.dat"
            new_file2_name = sim_name + "_restartBC.dat"
            #print(new_file1_name)

            # Define full paths for copied and renamed files
            new_file1_path = os.path.join(current_sim_folder, new_file1_name)
            new_file2_path = os.path.join(current_sim_folder, new_file2_name)

            # Copy the files
            shutil.copy(file1, new_file1_path)
            shutil.copy(file2, new_file2_path)

            print(Fore.GREEN + f"   [SUCCESS] air_5 restart files copied and renamed successfully")

        except Exception as e:
            print(Fore.RED + f"   [ERROR] The air_5 restart files couldn't be copied or renamed: {e}")

    else:
        print(Fore.RED + f"   [ERROR] air_5 restart files were not found in: {stagline_restart_simulations_path}")


    if (mixture in file1):
        print(Fore.WHITE + f"   [INFO] No columns needs to be added in restart file!")
    else:

        # Read the file1 into a DataFrame
        df = pd.read_csv(new_file1_path, delim_whitespace=True, header=None)

        # Define the new columns with the required value
        new_col1 = np.full((df.shape[0],), 1.0000000000E-10)
        new_col2 = np.full((df.shape[0],), 1.0000000000E-10)

        # Insert the new columns at the beginning
        df = pd.concat([pd.DataFrame({0: new_col1, 1: new_col2}), df], axis=1)

        # Save the modified file, preserving the format and without headers
        df.to_csv(new_file1_path, sep=" ", index=False, header=False, float_format="%.10E")

        # Read the file2 into a DataFrame
        df = pd.read_csv(new_file2_path, delim_whitespace=True, header=None)

        # Define the new columns with the required value
        new_col1 = np.full((df.shape[0],), 1.0000000000E-10)
        new_col2 = np.full((df.shape[0],), 1.0000000000E-10)

        # Insert the new columns at the beginning
        df = pd.concat([pd.DataFrame({0: new_col1, 1: new_col2}), df], axis=1)

        # Save the modified file, preserving the format and without headers
        df.to_csv(new_file2_path, sep=" ", index=False, header=False, float_format="%.10E")

def MPPSpeciesDensities(temperature,pressure,mixture):
    command = f"mppequil -T {temperature} -P {pressure} -s 3 {mixture}"
    result = subprocess.check_output(command, shell=True).decode("utf-8")
    
    # Extract the last line and split values while filtering out empty strings
    density_values = result.strip().split("\n")[-1].split()

    # Remove any empty strings from the list
    density_values = [val for val in density_values if val.strip()]
    
    return density_values

def MPPMassFractions(temperature,pressure,mixture):
    command = f"mppequil -T {temperature} -P {pressure} -s 2 {mixture}"
    result = subprocess.check_output(command, shell=True).decode("utf-8")
    
    # Extract the last line and split values while filtering out empty strings
    mass_fractions_values = result.strip().split("\n")[-1].split()

    # Remove any empty strings from the list
    mass_fractions_values = [val for val in mass_fractions_values if val.strip()]
    
    return mass_fractions_values

def MPPMixtureDensity(temperature,pressure,mixture):
    command = f"mppequil -T {temperature} -P {pressure} -m 3 {mixture}"
    result = subprocess.check_output(command, shell=True).decode("utf-8")
    
    # Extract the last line and split values while filtering out empty strings
    result = result.strip().split("\n")[-1].split()

    # Remove any empty strings from the list
    density_value = [val for val in result if val.strip()]

    density_value = [float(val) for val in density_value]
    
    return density_value[0]

def MPPViscosity(temperature,pressure,mixture):
    command = f"mppequil -T {temperature} -P {pressure} -m 32 {mixture}"
    result = subprocess.check_output(command, shell=True).decode("utf-8")
    
    # Extract the last line and split values while filtering out empty strings
    result = result.strip().split("\n")[-1].split()

    # Remove any empty strings from the list
    viscosity = [val for val in result if val.strip()]

    viscosity = [float(val) for val in viscosity]
    
    return viscosity[0]

def replace_inputs(lines, keyword, new_value):
        for i, line in enumerate(lines):
            if line.strip() == keyword:
                lines[i + 1] = f"{new_value}\n"
                break

def replace_init_cond(lines,init_cond):

        static_pressure = init_cond.get("Pressure")
        densities = init_cond.get("Densities")
        mass_fractions = init_cond.get("Mass fractions")
        velocities = init_cond.get("Velocities")
        temperatures = init_cond.get("Temperatures")

        for i, line in enumerate(lines):
            if line.strip() == "uin dv_dyin Tin pin yiin":
                    
                # Writing uin
                a = i + 1
                lines[a] = str(velocities[0]) + "\n"
                
                a += 1

                # Writing dv_dyin
                lines[a] = "0\n"

                a += 1

                # Writing Tin
                lines[a] = str(temperatures[0]) + "\n"

                a += 1

                # Writing static pressure 
                lines[a] = str(static_pressure) + "\n"

                a += 1

                # Writing mass fractions
                for value in mass_fractions:
                    lines[a] = str(value) + "\n"
                    a += 1  


            if line.strip() == "# Physical variable values (species densities, velocity components and temperatures) ":

                density_counter = 0

                for j, value in enumerate(densities):
                    
                    a = i + 1 + j

                    lines[a] = str(value) + "\n"

                    density_counter = j
                    

                velocity_counter = 0
                
                for j, value in enumerate(velocities):

                    a = i + 2 + density_counter + j

                    lines[a] = str(value) + "\n"

                    velocity_counter = j
                    
                
                for j, value in enumerate(temperatures):

                    a = i + 3 + density_counter + velocity_counter + j

                    lines[a] = str(value) + "\n"

def PdynToVelocity(pdyn, rho, mu, R):
    def equation(U):
        return pdyn - (0.5 * rho * U**2 * (1 + (6 / (((rho * U * R) / mu) + 0.455 * np.sqrt((rho * U * R) / mu)))))

    # Provide an initial guess for U
    U_initial_guess = np.sqrt(2 * pdyn / rho)

    # Solve using Newton's method
    solution_U = newton(equation, U_initial_guess)

    return solution_U

def NumberToFortranNotation(value):
    """Convert a decimal number to Fortran-style double precision scientific notation."""
    if value == 0:
        return "0.0d0"  # Special case for zero
    
    exponent = int(f"{value:.1e}".split('e')[1])  # Extract exponent
    mantissa = value / (10 ** exponent)  # Normalize mantissa
    
    return f"{mantissa:.1f}d{exponent}"

def RunStagline(sim_folder_path,sim_name,stagline_exe_path,CFL_range,Iter,cfl_inter,res_plot_visu):
    print(Fore.BLUE + "[STEP] Running Stagline simulation...")
    
    if not os.path.isdir(sim_folder_path):
        print(Fore.RED + f"[ERROR] The folder '{sim_folder_path}' does not exist!")
        return
    
    os.chdir(sim_folder_path)
    print(Fore.WHITE + f"   [INFO] Directory changed to: {os.getcwd()}")

    command = stagline_exe_path
    process = subprocess.Popen(command, shell=True)  # Run asynchronously

    if(res_plot_visu == True):
        plt.ion()  # Interactive mode ON
        fig, ax = plt.subplots()
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Residual")
        ax.set_title("Convergence Plot")
        line, = ax.plot([], [], 'b-', label='Convergence')
        ax.legend()

    convergence_file_name = sim_name + "_convergence.dat"
    cfl_file_name = "cfl"

    data_file = os.path.join(sim_folder_path, convergence_file_name)
    cfl_file = os.path.join(sim_folder_path, cfl_file_name)

    time.sleep(3)  # Wait for 2 seconds before starting to print

    processed_cfl = []
    
    while process.poll() is None:  # While the process is running

        if cfl_inter == ".TRUE.":
            try:
                # Check if the data file exists and is not empty
                if os.path.exists(data_file) and os.stat(data_file).st_size > 0:
                    with open(data_file, 'r') as f:
                        lines = [line for line in f.readlines() if not line.startswith('#')]

                    if lines:
                        # Load data while skipping the first two rows
                        file_convergence = np.loadtxt(data_file, skiprows=2)

                        if file_convergence.ndim > 0 and file_convergence.size > 0:
                            iterations = file_convergence[:, 0]

                            # Update CFL if a matching iteration is found
                            for i, iter_value in enumerate(Iter):
                                if iter_value in iterations and iter_value not in processed_cfl:
                                    new_cfl = NumberToFortranNotation(CFL_range[i])
                                    with open(cfl_file, 'w') as f:
                                        f.write(new_cfl)
                                        f.flush()  # Forces the write
                                    print(Fore.CYAN + f"[CFL UPDATE] Updated CFL to {new_cfl.strip()} at iteration {iter_value}")

                                    # Prevent redundant updates
                                    processed_cfl.append(Iter[i]) 
            except Exception as e:
                print(Fore.YELLOW + f"[WARNING] Error updating CFL: {e}")

        if res_plot_visu == True:
            # Update residuals plot
            try:
                if os.path.exists(data_file):
                    data = np.loadtxt(data_file)
                    if data.size > 0 and data.shape[1] >= 5:
                        x, y = data[:, 0], data[:, 4]  # Column 1 (x) and Column 5 (y)
                        line.set_xdata(x)
                        line.set_ydata(y)
                        ax.relim()
                        ax.autoscale_view()
                        plt.draw()
                        plt.pause(1)  # Refresh the plot every second
            except Exception as e:
                print(Fore.YELLOW + f"[WARNING] Could not update residual plot: {e}")

            time.sleep(1)  # Prevents excessive CPU usage



    process.wait()  # Wait for the simulation to finish

    processed_cfl.clear()

    if process.returncode == 0:
        print(Fore.GREEN + "[SUCCESS] Stagline simulation completed successfully.")
    else:
        print(Fore.RED + "[ERROR] Stagline simulation failed!")