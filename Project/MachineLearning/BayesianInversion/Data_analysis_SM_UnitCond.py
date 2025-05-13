import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az # type: ignore

from Library.Utils_BayesianInversion import *

# ========================================================================
 
# ------------------------------------
# | STAGLINE SURROGATE MODEL LOADING |
# ------------------------------------

model_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/KRGModel/SM_Gamma_Log"

sm_q = LoadModel(model_path)

# ========================================================================
 
# -----------------------------
# | EXPERIMENTAL DATA LOADING |
# -----------------------------

# Path to test campaign data
CSV_path = "/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Data/tests-mass-flow-Justin.xlsx"

# Columns to check for wrong data
columns_to_check = ["Pressure[mbar]", "Pitot[Pa]", "T [K] (x = 375mm, r = 0mm)", "HeatFlux(HS50mm)[kW/m2]"]

Pstat,massflow,power,HF,off_set_HF,Pdyn,off_set_Pdyn,temperature = CSVReader(CSV_path,columns_to_check)

# Target static pressure 
# ----------------------
target_pressure = [1500,5000,10000]
target_massflow = 16
tolerance_press = 300
tolerance_mf = 3

for targ_p in target_pressure:

    # Gathering data for bayesian inversion
    # -------------------------------------
    Pdyn_test = []
    Pstat_test = []
    T_test = []
    HF_test = []
    off_set_Pdyn_test = []
    off_set_HF_test = []



    for i,pstat in enumerate(Pstat):

        if pstat >= targ_p - tolerance_press and pstat <= targ_p + tolerance_press:
            
            if massflow[i] >= target_massflow - tolerance_mf and massflow[i] <= target_massflow + tolerance_mf:

                Pdyn_test.append(Pdyn[i])
                Pstat_test.append(Pstat[i])
                T_test.append(temperature[i])
                HF_test.append(HF[i])
                off_set_Pdyn_test.append(off_set_Pdyn[i])
                off_set_HF_test.append(off_set_HF[i])

    HF_test = np.array(HF_test)/1000

    # --------------------------------------------
    # | SURROGATE MODEL PREDICTIONS & PLOTTING   |
    # --------------------------------------------

    plt.figure(figsize=(8, 6))

    # Plot experimental data once before the loop
    plt.plot(
        T_test, HF_test,
        label="Experimental data",
        marker='D',
        linestyle='--',
        color='black'
    )

    for i in range(len(Pdyn_test)):
        # === Load the posterior samples from a .nc file ===
        trace_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Simulations/Gamma_priors/UniformG_UnitCond_Pstat={targ_p/100}/T={T_test[i]}/Posteriors/Thinned_trace.nc"  
        idata = az.from_netcdf(trace_path)

        # === Extract GN and GO samples ===
        GN_samples = idata.posterior["GN"].stack(sample=("chain", "draw")).values.flatten()
        GO_samples = idata.posterior["GO"].stack(sample=("chain", "draw")).values.flatten()

        # === Randomly select 1000 points ===
        n_points = 5000
        if len(GN_samples) < n_points:
            raise ValueError(f"Not enough samples: requested {n_points}, but got {len(GN_samples)}")

        indices = np.random.choice(len(GN_samples), size=n_points, replace=False)
        GN_selected = GN_samples[indices]
        GO_selected = GO_samples[indices]

        # === Predict for each sample ===
        predictions = []
        input_pairs = []
        for gn, go in zip(GN_selected, GO_selected):
            XV = np.array([Pdyn_test[i], Pstat_test[i], T_test[i], gn, go]).reshape(1, -1)
            prediction = sm_q.predict_values(XV)[0, 0] / 1000  # kW/m²
            predictions.append(prediction)
            input_pairs.append((gn, go))

        predictions = np.array(predictions)
        mean_pred = predictions.mean()
        std_pred = predictions.std() * 1.96

        # === Get the GN/GO pair closest to the mean prediction ===
        idx_closest = np.abs(predictions - mean_pred).argmin()
        gn_mean = 10**(input_pairs[idx_closest][0])
        go_mean = 10**(input_pairs[idx_closest][1])

        # === Plot prediction with error bar and annotate GN/GO ===
        plt.errorbar(
            T_test[i],
            mean_pred,
            yerr=std_pred,
            fmt='o',
            color='tab:blue',
            ecolor='lightblue',
            elinewidth=2,
            capsize=4,
            label="Predicted (mean ± 95% Conf. Int.)" if i == 0 else ""
        )

        plt.text(
            T_test[i] + 50, mean_pred + 200,
            f"GN={gn_mean:.5f}\nGO={go_mean:.5f}",
            fontsize=8,
            color="blue",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="blue", facecolor="white", alpha=0.7)
        )
    
    # Save plot
    output_path = f"Data_analysis/Unit_cond/Gamma_priors/UniformG_Pred_vs_Exp_Pstat={targ_p}.jpeg"
    plt.legend()
    plt.title(f"Experimental vs predicted HF - {targ_p/100} mbar")
    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()


        


            
