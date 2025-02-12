import shutil
import os
import subprocess
import numpy as np
import sympy as sp
from colorama import Fore, Style, init

init(autoreset=True)

def UserInputsProcessing(Txin, Tyin, pc, mixture, pdyn, uin, vin, R, n_species, cfl_val, Twall, cfl_inter, cfl_adaptive):
    # Initialisation of the inputs and init_cond dictionary
    # -----------------------------------------------------
    print(Fore.BLUE + "[STEP] Processing user inputs")
    
    inputs = {}
    init_cond = {}
    
    # Check for type of inputs (Pdyn or velocity)
    # -------------------------------------------
    if pdyn is not None and uin is None and vin is None:
        print(Fore.WHITE + "   [INFO] Calculating velocity from dynamic pressure")
        
        # Computing the mixture density
        rho = MPPMixtureDensity(Txin, pc, mixture)
        print(Fore.WHITE + f"   [INFO] Computed density: {rho}")

        # Computing the mixture dynamic viscosity
        mu = MPPViscosity(Txin, pc, mixture)
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
    sim_name = f"sim_Pc={pc}_T={Txin}_U={np.abs(uin):.1f}_mix={mixture}"
    
    # Processing the rest of the inputs
    # ---------------------------------
    inputs.update({
        "Simulation_Name": sim_name,
        "Mixture": mixture,
        "Number_of_species": n_species,
        "CFL_number": cfl_val,
        "Twall": Twall,
        "Inter_CFL": cfl_inter,
        "Adaptive_CFL": cfl_adaptive
    })
    
    # Gathering the species density
    density = MPPSpeciesDensities(Txin, pc, mixture)
    print(Fore.WHITE + f"   [INFO] Computed species densities: {density}")
    
    # Initial conditions values
    init_cond["Densities"] = density
    init_cond["Temperatures"] = np.array([Txin, Tyin])
    
    return inputs, init_cond, sim_name

def InputFileGenerator(stagline_simulations_path, input_template_path, sim_name, inputs, init_cond):
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
    mixture = inputs.get("Mixture")
    
    try:
        shutil.copy(os.path.join(input_template_path, f"Example_input_{mixture}"), sim_folder_path)
        shutil.copy(os.path.join(input_template_path, "mesh.dat"), sim_folder_path)
        print(Fore.WHITE + "   [INFO] Reference files successfully copied.")
    except FileNotFoundError:
        print(Fore.RED + "[ERROR] One or more reference files are missing!")
        return None
    
    # Renaming the input file
    sim_input = os.path.join(sim_folder_path, f"Example_input_{mixture}")
    renamed_input = os.path.join(sim_folder_path, "input")
    os.rename(sim_input, renamed_input)
    
    print(Fore.WHITE + "   [INFO] Modification of the input file with user inputs")

    with open(renamed_input, "r") as input_file:
        lines = input_file.readlines()
    
    for key, value in inputs.items():
        replace_inputs(lines, key, str(value))
    
    replace_init_cond(lines, init_cond)
    
    with open(renamed_input, 'w') as input_file:
        input_file.writelines(lines)
    
    print(Fore.GREEN + "   [SUCCESS] Input file successfully modified.")
    return sim_folder_path

def MPPSpeciesDensities(temperature,pressure,mixture):
    command = f"mppequil -T {temperature} -P {pressure} -s 3 {mixture}"
    result = subprocess.check_output(command, shell=True).decode("utf-8")
    
    # Extract the last line and split values while filtering out empty strings
    density_values = result.strip().split("\n")[-1].split()

    # Remove any empty strings from the list
    density_values = [val for val in density_values if val.strip()]
    
    return density_values

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

        densities = init_cond.get("Densities")
        velocities = init_cond.get("Velocities")
        temperatures = init_cond.get("Temperatures")

        for i, line in enumerate(lines):
            if line.strip() == "rhoiin uin vin Tin Tvin":

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

            if line.strip() == "# Physical variable values (species densities, velocity components and temperatures)":

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
    U = sp.Symbol('U', real=True, positive=True)
    equation = sp.Eq(pdyn, (1/2) * rho * U**2 * (1 + (6 / (((rho * U * R) / mu) + 0.455 * sp.sqrt((rho * U * R) / mu)))))
    solution_U = sp.solve(equation, U)

    return solution_U[0]

def RunStagline(sim_folder_path):
    # Change directory to run stagline
    print(Fore.BLUE + "[STEP] Running Stagline simulation...")
    
    if not os.path.isdir(sim_folder_path):
        print(Fore.RED + f"[ERROR] The folder '{sim_folder_path}' does not exist!")
        return
    
    os.chdir(sim_folder_path)
    print(Fore.WHITE + f"   [INFO] Directory changed to: {os.getcwd()}")
    
    command = "../../bin/stagline"
    result = subprocess.run(command, shell=True)
    
    if result.returncode == 0:
        print(Fore.GREEN + "[SUCCESS] Stagline simulation completed successfully.")
    else:
        print(Fore.RED + "[ERROR] Stagline simulation failed!")
