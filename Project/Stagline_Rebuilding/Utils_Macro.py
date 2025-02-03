import shutil
import os
import subprocess
import numpy as np

def InputFileGenerator(sim_name, inputs, init_cond):
    """
    Generate and modify the input file for a Stagline simulation.
    """
    print("\033[94mStarting input file modification for Stagline simulation...\033[0m")
    
    folder_path = "/home/jpe/VKI/Project/Stagline/Simulations"
    
    if os.path.isdir(folder_path):
        print(f"[INFO] Stagline directory found: {folder_path}")
    else:
        print("\033[91m[ERROR] Simulation directory not found!")
        return None
    
    # Creating simulation folder
    sim_folder_path = os.path.join(folder_path, sim_name)
    if not os.path.isdir(sim_folder_path):
        os.makedirs(sim_folder_path)
        print(f"[INFO] Created simulation folder: {sim_folder_path}")
    else:
        print(f"\033[93m[WARNING] Simulation folder already exists: {sim_folder_path}\033[0m")
    
    # Copying necessary files from macro folder
    mixture = inputs.get("Mixture")
    macro_folder = "/home/jpe/VKI/Project/Stagline/Macro_folders" 
    
    try:
        shutil.copy(os.path.join(macro_folder, f"Example_input_{mixture}"), sim_folder_path)
        shutil.copy(os.path.join(macro_folder, "mesh.dat"), sim_folder_path)
        print("[INFO] Reference files successfully copied.")
    except FileNotFoundError:
        print("\033[91m[ERROR] One or more reference files files are missing!\033[0m")
        return None
    
    # Modifying the input file
    sim_input = os.path.join(sim_folder_path, f"Example_input_{mixture}")
    renamed_input = os.path.join(sim_folder_path, "input")
    os.rename(sim_input, renamed_input)
    
    with open(renamed_input, "r") as input_file:
        lines = input_file.readlines()
    
    for key, value in inputs.items():
        replace_inputs(lines, key, str(value))
    
    replace_init_cond(lines, init_cond)
    
    with open(renamed_input, 'w') as input_file:
        input_file.writelines(lines)
    
    print("\033[92m[SUCCESS] Input file successfully modified.\033[0m\n")
    return sim_folder_path

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
    
def MPPDensities(temperature,pressure,mixture):
    command = f"mppequil -T {temperature} -P {pressure} -s 3 {mixture}"
    result = subprocess.check_output(command, shell=True).decode("utf-8")
    
    # Extract the last line and split values while filtering out empty strings
    density_values = result.strip().split("\n")[-1].split()

    # Remove any empty strings from the list
    density_values = [val for val in density_values if val.strip()]
    
    return density_values

def RunStagline(sim_folder_path):

    # Change directory to run stagline
    if os.path.isdir(sim_folder_path):
        os.chdir(sim_folder_path)
        print(f"    Directory changed to: {os.getcwd()}")
    else:
        print("\033[91m" + f"Error: The folder '{sim_folder_path}' does not exist!" + "\033[0m")

    command = "../../bin/stagline"
    result = subprocess.run(command, shell=True)
    
        