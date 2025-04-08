import numpy as np
from colorama import Fore, Style, init
from smt.surrogate_models import KRG
import sys
import pandas as pd
import matplotlib.pylab as plt
from matplotlib.lines import Line2D


# Initialize colorama for colored terminal output
init(autoreset=True)

def LoadModel(model_path):

    print(Fore.BLUE + "[STEP] Loading the surrogate model")

    print(Fore.WHITE + f"---> [INFO] Loading Model: '{model_path}' ...")
    try:
        sm_q = KRG.load(model_path)
        sm_q.options["print_global"] = False
        print(Fore.GREEN + f"---> [SUCCESS] Model '{model_path}' loaded successfully!")
        return sm_q
    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Failed to load model '{model_path}': {e}")
        sys.exit(1)

def LoadDataEnrico(StagInput_path,StagRes_path):

    print(Fore.BLUE + "[STEP] Loading data from csv files")

    print(Fore.WHITE + f"---> [INFO] Loading Stagline inputs ...")
    # Loading Stagline inputs from csv file
    df_Inputs = pd.read_excel(StagInput_path)

    print(Fore.WHITE + f"---> [INFO] Loading Stagline data ...")
    # Loading Stagline results from csv file
    df_Res = pd.read_excel(StagRes_path)

    # Loop on Stagline results
    full_data = []
    for index, row in df_Res.iterrows():

        # Extract results
        Pc     = row["Pc"]
        T      = row["T"]
        Qstag  = row["Q_stag"]
        gamma  = row["gamma"]

        # Find matching lines in input csv
        line_match = df_Inputs[(df_Inputs["Pc"] == Pc) & (df_Inputs["T"] == T)]

        # Extract inputs values
        if not line_match.empty:
            line = line_match.iloc[0]  

            Pdyn = line["Pdyn"]

            # Adding data into the list
            full_data.append({
                "Pc": Pc,
                "T": T,
                "Pdyn": Pdyn,
                "gamma": gamma,
                "Q_stag": Qstag
            })
        else:
            print(Fore.RED + f"---> [ERROR] No matching Pdyn for Pc = {Pc} and T = {T}")
            sys.exit(1)

    # Creating the full dataframe

    df_FullData = pd.DataFrame(full_data)

    print(Fore.GREEN + f"---> [SUCCESS] Validation data successfully extracted !")

    return df_FullData
        
def RunSM(sm_q,df_FullData):

    print(Fore.BLUE + "[STEP] Model predictions")

    print(Fore.WHITE + "---> [INFO] Running predictions on validation data ...")

    Results = []

    for index, row in df_FullData.iterrows():

        # Gathering Data to run SM
        Pc = row["Pc"]
        T = row["T"]
        Pdyn = row["Pdyn"]
        gamma = row["gamma"]
        Q_stag = row["Q_stag"]

        # Gamma must be log !!
        if gamma == 0:
            gamma_corr = 1e-4
            gamma_log = np.log10(gamma_corr)
        else:
            gamma_log = np.log10(gamma)

        # Preparing the Data vector for the SM (shape = (5,1))
        XV = np.array([Pdyn,Pc*100,T,gamma_log,gamma_log]).reshape((1,5))
        try:
            # Y predictions
            YV_K= sm_q.predict_values(XV)/1000
            YV_K = YV_K.item()
            # Get prediction uncertainty (variance)
            #YV_K_var = sm_q.predict_variances(XV) / 1000000  # Convert from W² to kW²
            #YV_K_std = np.sqrt(YV_K_var)  # Convert variance to standard deviation
        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Prediction failed: {e}")
            sys.exit(1)

        # Computing error
        err = np.abs(Q_stag - YV_K)

        # Adding data into the list
        Results.append({
            "Pdyn": Pdyn,
            "Pc": Pc,
            "T": T,
            "gamma": gamma,
            "Q_stag": Q_stag,
            "Q_pred": YV_K,
            "Error": err
        })
        
    print(Fore.GREEN + "---> [SUCCESS] Model Predictions successful!")

    # Creating Dataframe with all the results
    print(Fore.WHITE + "---> [INFO] Creating Dataframe with all the results ...")
    df_Results = pd.DataFrame(Results)

    # Removing values for gamma = 1
    df_nrmse_gamma_1 = df_Results[df_Results["gamma"] != 1]
    df_nrmse = df_Results

    # Compute NRMSE without gamma = 1
    err_eval = 0
    err_eval = np.sum((df_nrmse["Q_stag"] - df_nrmse["Q_pred"])**2)

    nrmse = np.sqrt(err_eval / len(df_nrmse)) * 100 / (df_nrmse["Q_pred"].max() - df_nrmse["Q_pred"].min())

    # Ensure `rmse` is a float before printing
    print(Fore.YELLOW + f"---> [RESULT] NRMSE Error: {float(nrmse):.4f}%")

    # Compute NRMSE with gamma = 1
    err_eval = 0
    err_eval = np.sum((df_nrmse_gamma_1["Q_stag"] - df_nrmse_gamma_1["Q_pred"])**2)

    nrmse = np.sqrt(err_eval / len(df_nrmse_gamma_1)) * 100 / (df_nrmse_gamma_1["Q_pred"].max() - df_nrmse_gamma_1["Q_pred"].min())

    # Ensure `rmse` is a float before printing
    print(Fore.YELLOW + f"---> [RESULT] NRMSE Error: {float(nrmse):.4f}% (without gamma = 1)")

    try:
        path_restults = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Validation_Enrico/Results/Validation_Results.xlsx"
        df_Results.to_excel(path_restults, index=False)
        print(Fore.GREEN + f"---> [SUCCESS] Results successfully dumped in {path_restults}!")

    except Exception as e:
        print(Fore.RED + f"---> [ERROR] Results dumping failed: {e}")
    

    return df_Results

def Plot(Pc,df_Results):
    
    # Gather results based on Pc
    df_Pc_sorted = df_Results[df_Results["Pc"] == Pc]

    # Get different values of gamma 
    gamma_val = df_Pc_sorted["gamma"].unique()

    plt.figure()

    print(Fore.WHITE + f"---> [INFO] Creating the figures for PC = {Pc}mbar ...")
    for gamma in gamma_val:

        if gamma == 1:

            plt.figure()

            # Sorting dataframe for gamma
            df_gamma_sorted = df_Pc_sorted[df_Pc_sorted["gamma"] == gamma]

            # Get Results
            T = df_gamma_sorted["T"]
            Q_pred = df_gamma_sorted["Q_pred"]
            Q_stag = df_gamma_sorted["Q_stag"]

            plt.plot(T, Q_pred, linestyle='--', marker='o', color='blue')
            plt.plot(T, Q_stag, linestyle='-',  marker='o', color='black')
            #plt.text(T.values[-1], Q_pred.values[-1], fr'$\gamma = {gamma}$', fontsize=12, rotation=45)
            plt.text(T.values[-1], Q_stag.values[-1], fr'$\gamma = {gamma}$', fontsize=12, rotation=45)


            plt.text(0.02, 0.90, '--  = Prediction', transform=plt.gca().transAxes, fontsize=11, color = 'blue')
            plt.text(0.02, 0.95, '—   = Stagline',   transform=plt.gca().transAxes, fontsize=11, color = 'black')

            # Titres et axes si tu veux
            plt.xlabel("Temperature [K]")
            plt.ylabel("Heat flux [kW]")
            plt.grid(True)
            plt.title(f"Results for Pc = {Pc} mbar")
            plt.savefig(f"Results/Validation_Pc={Pc}_gamma={gamma}.jpeg",format="jpeg",dpi=300, bbox_inches='tight')
            plt.close()
        else:

            # Sorting dataframe for gamma
            df_gamma_sorted = df_Pc_sorted[df_Pc_sorted["gamma"] == gamma]

            # Get Results
            T = df_gamma_sorted["T"]
            Q_pred = df_gamma_sorted["Q_pred"]
            Q_stag = df_gamma_sorted["Q_stag"]

            plt.plot(T, Q_pred, linestyle='--', marker='o', color='blue')
            plt.plot(T, Q_stag, linestyle='-',  marker='o', color='black')
            #plt.text(T.values[-1], Q_pred.values[-1], fr'$\gamma = {gamma}$', fontsize=12, rotation=45)
            plt.text(T.values[-1], Q_stag.values[-1], fr'$\gamma = {gamma}$', fontsize=12, rotation=45)


    plt.text(0.02, 0.90, '--  = Prediction', transform=plt.gca().transAxes, fontsize=11, color = 'blue')
    plt.text(0.02, 0.95, '—   = Stagline',   transform=plt.gca().transAxes, fontsize=11, color = 'black')

    # Titres et axes si tu veux
    plt.xlabel("Temperature [K]")
    plt.ylabel("Heat flux [kW]")
    plt.grid(True)
    plt.title(f"Results for Pc = {Pc} mbar")
    plt.savefig(f"Results/Validation_Pc={Pc}.jpeg",format="jpeg",dpi=300, bbox_inches='tight')
    plt.close()

def PerformanceEval(df_results):

    df_results = df_results[df_results["gamma"] != 1]

    print(Fore.WHITE + f"---> [INFO] Plotting performance evaluation graph ...")

    YV = df_results["Q_stag"]
    YV_K = df_results["Q_pred"]
    # Create the plot
    plt.figure(figsize=(6, 6))
    # Plot the confidence interval with smooth shading
    plt.plot([np.min(YV), np.max(YV)], [np.min(YV), np.max(YV)], 'r--', label="Perfect Fit")
    #plt.fill_between(
    #    Y_sorted, YV_K_lower, YV_K_upper,
    #    color="orange", alpha=0.4, label=r"95% Confidence Interval $[kW^2]$"
    #)
    # Scatter plot for actual and predicted values
    plt.scatter(YV, YV_K, label="Stagline data vs predictions", alpha=0.7)
    # Add labels and title
    plt.legend(loc="upper left")
    plt.xlabel(r"$q_{wall}$ from Stagline [kW]", fontsize=14)
    plt.ylabel(r"$q_{wall}$ predicted [kW]", fontsize=14)
    plt.title("Model Prediction Performance")
    # Save figure in high resolution
    plt.savefig("Results/ModelPerformance.jpeg", format='jpeg', dpi=300, bbox_inches='tight')





        









    