  
from Utils_StaglineSupFastRun import *


def StaglineFastRunSupersonic(Tin,Tvin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart,stagline_simulations_path,input_template_path):
    # -----------------------------------------------
    # | START OF THE SCRIPT FOR STAGLINE SIMULATION |
    # -----------------------------------------------
    
    # Processing user inputs
    # ----------------------
    inputs,init_cond,sim_name,mixture = UserInputsProcessing(Tin,Tvin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart)
    
    # Runing script for input file modification
    # -----------------------------------------
    sim_folder_path = InputFileGenerator(stagline_simulations_path,input_template_path,sim_name,inputs,init_cond,mixture)
    
    # Stagline simulation
    # -------------------
    RunStagline(sim_folder_path,sim_name)