from Utils_StaglineSubFastRun import *

def StaglineFastRunSubsonic(Tin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart,stagline_simulations_path,input_template_path,stagline_exe_path,CFL_range,Iter,air_5_restart):
    # -----------------------------------------------
    # | START OF THE SCRIPT FOR STAGLINE SIMULATION |
    # -----------------------------------------------
    
    # Processing user inputs
    # ----------------------
    inputs,init_cond,sim_name,mixture = UserInputsProcessing(Tin,pc,mixture,pdyn,uin,vin,R,n_species,cfl_val,Twall,cfl_inter,cfl_adaptive,Log_CFL,residual,restart)
    
    # Runing script for input file modification
    # -----------------------------------------
    sim_folder_path = InputFileGenerator(stagline_simulations_path,input_template_path,sim_name,inputs,init_cond,mixture,air_5_restart)
    
    # Stagline simulation
    # -------------------
    RunStagline(sim_folder_path,sim_name,stagline_exe_path,CFL_range,Iter,cfl_inter)