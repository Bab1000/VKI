import shutil
import os
import subprocess
import numpy as np

def InputFileGenerator(sim_name, inputs,init_cond):
    # Entering stagline folder

    print("\033[94mModification of the input file for stagline simulation\033[0m")

    folder_path = "/home/jpe/VKI/Project/Stagline/Simulations"

    if os.path.isdir(folder_path):
        print(f"    The folder {folder_path} exists ...")
    else:
       print("\033[91m" + f"The folder {folder_path} does not exist !" + "\033[0m")


    
    # Creating working folder
    sim_folder = str(sim_name)
    sim_folder_path = folder_path + '/' + sim_folder 

    if not os.path.isdir(sim_folder_path):
        os.makedirs(sim_folder_path)
        print(f"    The folder {sim_folder_path} has been created ...")
    else:
        print(f"The folder {sim_folder_path} already exists ...")

    # Cpying simulation files from the macro folder
    mixture = inputs.get("Mixture")
    macro_folder = "/home/jpe/VKI/Project/Stagline/Macro_folders" 
    macro_input = macro_folder + "/Example_input_" + mixture
    macro_mesh = macro_folder + "/mesh.dat"
    shutil.copy(macro_input,sim_folder_path)
    shutil.copy(macro_mesh,sim_folder_path)

    # Modification to the input file
    sim_input = sim_folder_path + "/Example_input_" + mixture
    os.rename(sim_input,sim_folder_path + "/input")
    sim_input = sim_folder_path + "/input"

    with open(sim_input, "r") as input_file:
        lines = input_file.readlines()
                    
    for key, value in inputs.items():

        replace_inputs(lines, key, str(value))
    
    replace_init_cond(lines,init_cond)

    with open(sim_input, 'w') as input_file:
        input_file.writelines(lines)

    print("\033[92mThe input file has been successfully modified\033[0m")
    print("")

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
    
        