import matplotlib.pyplot as plt
import numpy as np
import arviz as az
import os
from colorama import Fore
from smt.surrogate_models import KRG
import sys
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

def plot_posteriors_individual(idata_dict, pressure_labels,
                               prior_range=(-4, 0),
                               output_path=".", filename_base="gamma_posteriors"):
    """
    Plots and saves individual posterior distribution plots (GN, GO, Gglob) per pressure.
    
    Parameters:
    - idata_dict: dict {pressure: (idata_ssg, idata_ug)}
    - pressure_labels: list of pressures to plot (e.g. [1500, 5000, 10000])
    - prior_range: x-range for all plots
    - output_path: folder where to save the figure
    - filename_base: base name for each output file (pressure will be appended)
    """
    x_grid = np.linspace(prior_range[0], prior_range[1], 500)
    prior_height = 1 / (prior_range[1] - prior_range[0])

    os.makedirs(output_path, exist_ok=True)

    for p in pressure_labels:
        idata_ssg, idata_ug = idata_dict[p]
        samples_ssg = idata_ssg.posterior.stack(sample=("chain", "draw"))
        samples_ug = idata_ug.posterior.stack(sample=("chain", "draw"))

        GN = samples_ssg["GN"].values.flatten()
        GO = samples_ssg["GO"].values.flatten()
        Gglob = samples_ug["Gglob"].values.flatten()

        kde_gn = gaussian_kde(GN)(x_grid)
        kde_go = gaussian_kde(GO)(x_grid)
        kde_gglob = gaussian_kde(Gglob)(x_grid)

        fig, ax = plt.subplots(figsize=(8, 5))

        # GN
        ax.plot(x_grid, kde_gn, color="#4A90E2", label=r"$\log_{10}(\gamma_N)$", linewidth=2)
        ax.fill_between(x_grid, kde_gn, color="#4A90E2", alpha=0.3)

        # GO
        ax.plot(x_grid, kde_go, color="#2ECC71", label=r"$\log_{10}(\gamma_O)$", linewidth=2)
        ax.fill_between(x_grid, kde_go, color="#2ECC71", alpha=0.3)

        # Gglob
        ax.plot(x_grid, kde_gglob, color="#9B59B6", label=r"$\log_{10}(\gamma_{glob})$", linewidth=2)
        ax.fill_between(x_grid, kde_gglob, color="#9B59B6", alpha=0.3)

        # Prior
        ax.axhline(prior_height, color="red", linestyle="--", linewidth=2, label="Uniform Prior")

        # Styling
        ax.set_xlabel(r"$\log_{10}(\gamma)$", fontsize=18)
        ax.set_ylabel("Probability Density", fontsize=18)
        ax.set_title(f"{p/100} mbar", fontsize=18)
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)
        ax.tick_params(labelsize=16)
        ax.legend(fontsize=14, loc='upper left')

        # Save
        filename = f"{filename_base}_{p/100}mbar.png"
        fig_path = os.path.join(output_path, filename)
        plt.tight_layout()
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"[SUCCESS] Saved posterior plot for {p/100} mbar → {fig_path}")

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

def MAP_plots_SSG(sm_model, samples, targ_p, Pstat_test, Pdyn_test, T_test,
              T_exp_list, HF_exp_list, HF_exp_err_list,
              output_path, global_output_path):
    """
    Loops over all posterior samples and plots MAP predictions vs experimental data.
    """

    print(Fore.BLUE + f"[STEP] Computing predictions for Pstat={targ_p/100} mbar")

    plt.figure(figsize=(8, 6))
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # === Plot experimental data ===
    plt.errorbar(
        T_exp_list, HF_exp_list,
        yerr=HF_exp_err_list,
        fmt='o',
        color='black',
        ecolor='gray',
        elinewidth=1.5,
        capsize=4,
        markersize=6,
        label="Experimental (± error)"
    )

    # === Setup pression statique ===
    pstat_color_dict = {1500: 'blue', 5000: 'green', 10000: 'red'}
    pstat_ecolor_dict = {1500: 'lightblue', 5000: 'lightgreen', 10000: 'lightcoral'}
    pstat_marker_dict = {1500: 's', 5000: 'D', 10000: '^'}  # carré, losange, triangle

    labels_done = set()

    for i in range(len(T_exp_list)):
        try:
            print(Fore.WHITE + f"---> [INFO] Computing MAP HF : Pstat = {targ_p/100} | Pdyn = {Pdyn_test[i]:.1f} | T = {T_test[i]:.1f} ...")

            try:
                Pdyn_samples = Pdyn_test
                Pstat_samples = Pstat_test
                T_samples = T_test
                G_samples = samples["Gglob"].values.flatten()
            except KeyError as e:
                print(Fore.RED + f"[ERROR] Missing variable in trace: {e}")
                continue

            n_points = 10000
            min_len = min(len(Pdyn_samples), len(Pstat_samples), len(T_samples), len(G_samples))
            if min_len < n_points:
                raise ValueError(f"Not enough samples: requested {n_points}, but got {min_len}")

            predictions = []
            for idx in range(min_len):
                XV = np.array([
                    Pdyn_samples[idx],
                    Pstat_samples[idx],
                    T_samples[idx],
                    G_samples[idx],
                    G_samples[idx]
                ]).reshape(1, -1)

                pred = sm_model.predict_values(XV)[0, 0] / 1000
                predictions.append(pred)

            predictions = np.array(predictions)
            map_hf = get_map(predictions, HF_exp_list[i], global_output_path, T_test[i], targ_p)
            std_pred = predictions.std() * 1.96

            # Récupérer la valeur moyenne de Pstat arrondie à la valeur connue
            avg_pstat = int(round(Pstat_samples.mean() / 100) * 100)  # ex: 1482 → 1500
            closest_pstat = min(pstat_color_dict.keys(), key=lambda x: abs(x - avg_pstat))

            label = f"MAP ± range (Pstat = {closest_pstat} Pa)"
            if label in labels_done:
                label = None
            else:
                labels_done.add(label)

            plt.errorbar(
                T_test[i],
                map_hf,
                yerr=std_pred,
                fmt=pstat_marker_dict[closest_pstat],
                color=pstat_color_dict[closest_pstat],
                ecolor=pstat_ecolor_dict[closest_pstat],
                elinewidth=2,
                capsize=4,
                markersize=7,
                alpha=0.9,
                label=label
            )

        except Exception as e:
            print(Fore.RED + f"[ERROR] MAP computation failed for point {i}: {e}")

    plt.legend()
    plt.title("Experimental vs MAP Heat Flux")
    plt.xlabel("T [K]")
    plt.ylabel("HF [kW/m2]")
    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

def get_map(predictions,HF_exp_val,global_output_path,T,targ_p):
    """
    Plots the KDE of the predicted HF values and highlights the MAP value.

    Parameters:
        predictions: array-like, predicted HF values (e.g., in kW/m²)
    """
    predictions = np.asarray(predictions)

    # Compute KDE and MAP
    kde = gaussian_kde(predictions)
    x_vals = np.linspace(predictions.min(), predictions.max(), 1000)
    y_vals = kde(x_vals)

    # Compute MAP (maximum a posteriori estimate)
    x_fine = np.linspace(predictions.min(), predictions.max(), 10000)
    map_val = x_fine[np.argmax(kde(x_fine))]

    return map_val

def MAP_comparison_plot_individual(sm_model,
                                   Pdyn_test, Pstat_test, T_test,
                                   HF_exp_list, HF_exp_err_list, global_output_path,
                                   filename_base="MAP_plot_SSG_vs_UG"):
    """
    Generates one figure per pressure bin (1500, 5000, 10000 Pa), showing:
    - Experimental data (black circles)
    - MAP prediction SSG (blue squares)
    - MAP prediction UG (green triangles)

    Each figure is saved separately with suffix _<pressure>mbar.
    """

    print(Fore.BLUE + "[STEP] MAP Comparison Plot — one figure per pressure")

    bins = [1500, 5000, 10000]
    colors = {'exp': 'black', 'ssg': '#1f77b4', 'ug': '#2ca02c'}
    markers = {'exp': 'o', 'ssg': 's', 'ug': '^'}
    alpha = 0.9
    cache = {}

    for bin_p in bins:
        fig, ax = plt.subplots()

        for i in range(len(T_test)):
            try:
                pstat_i = Pstat_test[i]
                bin_i = min(bins, key=lambda x: abs(x - pstat_i))
                if bin_i != bin_p:
                    continue  # skip points from other bins

                print(Fore.WHITE + f"---> [INFO] Bin {bin_p} Pa | Point {i}: T = {T_test[i]:.1f} K | Pstat = {pstat_i:.1f} Pa")

                # Load posterior samples if needed
                if bin_p not in cache:
                    file_ssg = os.path.join(global_output_path, f"Gamma_models/PBC_SSG_{int(bin_p/100)}mbar.nc")
                    file_ug  = os.path.join(global_output_path, f"Gamma_models/PBC_UG_{int(bin_p/100)}mbar.nc")
                    cache[bin_p] = {
                        "ssg": az.from_netcdf(file_ssg).posterior.stack(sample=("chain", "draw")),
                        "ug":  az.from_netcdf(file_ug).posterior.stack(sample=("chain", "draw"))
                    }

                samples_ssg = cache[bin_p]["ssg"]
                samples_ug = cache[bin_p]["ug"]

                # Plot experimental data
                ax.errorbar(
                    T_test[i], HF_exp_list[i],
                    yerr=HF_exp_err_list[i],
                    fmt=markers["exp"],
                    color=colors["exp"],
                    ecolor='gray',
                    elinewidth=1,
                    capsize=4,
                    markersize=6,
                    alpha=alpha
                )

                # SSG prediction
                preds_ssg = [
                    sm_model.predict_values(np.array([
                        Pdyn_test[i], Pstat_test[i], T_test[i],
                        samples_ssg["GN"][j], samples_ssg["GO"][j]
                    ]).reshape(1, -1))[0, 0] / 1000
                    for j in range(len(samples_ssg["GN"]))
                ]
                map_ssg = np.median(preds_ssg)
                std_ssg = np.std(preds_ssg) * 1.96

                ax.errorbar(
                    T_test[i], map_ssg,
                    yerr=std_ssg,
                    fmt=markers["ssg"],
                    color=colors["ssg"],
                    ecolor=colors["ssg"],
                    elinewidth=2,
                    capsize=4,
                    markersize=8,
                    alpha=alpha
                )

                # UG prediction
                preds_ug = [
                    sm_model.predict_values(np.array([
                        Pdyn_test[i], Pstat_test[i], T_test[i],
                        samples_ug["Gglob"][j], samples_ug["Gglob"][j]
                    ]).reshape(1, -1))[0, 0] / 1000
                    for j in range(len(samples_ug["Gglob"]))
                ]
                map_ug = np.median(preds_ug)
                std_ug = np.std(preds_ug) * 1.96

                ax.errorbar(
                    T_test[i], map_ug,
                    yerr=std_ug,
                    fmt=markers["ug"],
                    color=colors["ug"],
                    ecolor=colors["ug"],
                    elinewidth=2,
                    capsize=4,
                    markersize=8,
                    alpha=alpha
                )

            except Exception as e:
                print(Fore.RED + f"[ERROR] MAP computation failed at point {i}: {e}")

        # Styling
        ax.set_xlabel("T [K]", fontsize=16)
        ax.set_ylabel("HF [kW/m²]", fontsize=16)
        ax.set_title(f"{bin_p/100} mbar", fontsize=16)
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

        # Legend
        handles = [
            Line2D([0], [0], marker=markers['exp'], color=colors['exp'], linestyle='None', label='Experimental values', markersize=7),
            Line2D([0], [0], marker=markers['ssg'], color=colors['ssg'], linestyle='None', label='Predictions SSG', markersize=7),
            Line2D([0], [0], marker=markers['ug'],  color=colors['ug'],  linestyle='None', label='Predictions UG', markersize=7),
        ]
        ax.legend(handles=handles, fontsize=12, loc='upper left')

        # Save
        fig_name = f"{filename_base}_{bin_p}mbar.png"
        fig_path = os.path.join(global_output_path, fig_name)
        plt.tight_layout()
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(Fore.GREEN + f"[SUCCESS] Figure saved to {fig_path}")

def evaluate_model_comparison(sm_model,
                                           Pdyn_list, Pstat_list, T_list,
                                           HF_exp_list,
                                           global_output_path):
    """
    Compare MAP predictions from SSG and UG gamma models vs experimental heat flux.
    Loads .nc files based on each Pstat and prints intermediate steps.

    Returns:
    - y_pred_ssg, y_pred_ug : MAP predictions
    - error_dict            : dict with MAPE and MAE for both models
    """

    y_true = np.array(HF_exp_list)
    y_pred_ssg = []
    y_pred_ug = []

    cache = {}
    print(Fore.BLUE + "\n[STEP] Starting MAP prediction with dynamic .nc loading")
    print(f"[INFO] Number of test points: {len(HF_exp_list)}\n")

    for i in range(len(HF_exp_list)):
        print(Fore.CYAN + f"\n--- Point {i} ---")

        pstat_i = Pstat_list[i]
        pdyn_i = Pdyn_list[i]
        temp_i = T_list[i]
        hf_exp = HF_exp_list[i]

        print(Fore.WHITE + f"[INFO] Pstat = {pstat_i:.1f}, Pdyn = {pdyn_i:.1f}, T = {temp_i:.1f}, HF_exp = {hf_exp:.3f}")

        bin_pstat = int(round(pstat_i / 100) * 100)
        closest_bin = min([1500, 5000, 10000], key=lambda x: abs(x - bin_pstat))
        suffix = f"{int(closest_bin / 100)}mbar"

        print(Fore.YELLOW + f"[INFO] → Assigned to .nc file bin: {suffix}")

        if suffix not in cache:
            print(Fore.MAGENTA + f"[LOAD] Loading SSG and UG posteriors for {suffix}...")
            file_ssg = os.path.join(global_output_path, f"Gamma_models/PBC_SSG_{suffix}.nc")
            file_ug  = os.path.join(global_output_path, f"Gamma_models/PBC_UG_{suffix}.nc")

            samples_ssg = az.from_netcdf(file_ssg).posterior.stack(sample=("chain", "draw"))
            samples_ug  = az.from_netcdf(file_ug).posterior.stack(sample=("chain", "draw"))

            cache[suffix] = {"ssg": samples_ssg, "ug": samples_ug}
        else:
            print(Fore.GREEN + f"[CACHE] Using cached posterior for {suffix}")

        samples_ssg = cache[suffix]["ssg"]
        samples_ug = cache[suffix]["ug"]

        # === MAP prediction - SSG
        print(Fore.WHITE + "[COMPUTE] Predicting with SSG model...")
        preds_ssg = [
            sm_model.predict_values(np.array([
                pdyn_i, pstat_i, temp_i,
                samples_ssg["GN"][j], samples_ssg["GO"][j]
            ]).reshape(1, -1))[0, 0] / 1000
            for j in range(len(samples_ssg["GN"]))
        ]
        map_ssg = np.median(preds_ssg)
        y_pred_ssg.append(map_ssg)

        # === MAP prediction - UG
        print(Fore.WHITE + "[COMPUTE] Predicting with UG model...")
        preds_ug = [
            sm_model.predict_values(np.array([
                pdyn_i, pstat_i, temp_i,
                samples_ug["Gglob"][j], samples_ug["Gglob"][j]
            ]).reshape(1, -1))[0, 0] / 1000
            for j in range(len(samples_ug["Gglob"]))
        ]
        map_ug = np.median(preds_ug)
        y_pred_ug.append(map_ug)

        print(Fore.GREEN + f"[RESULT] MAP SSG = {map_ssg:.3f} | MAP UG = {map_ug:.3f}")

    # === Error computation
    print(Fore.BLUE + "\n[STEP] Final error metrics computation")

    y_pred_ssg = np.array(y_pred_ssg)
    y_pred_ug = np.array(y_pred_ug)

    def mape(y_true, y_pred):
        return np.mean(np.abs((y_pred - y_true) / np.abs(y_true))) * 100

    def mae(y_true, y_pred):
        return np.mean(np.abs(y_pred - y_true))

    mape_ssg = mape(y_true, y_pred_ssg)
    mae_ssg = mae(y_true, y_pred_ssg)

    mape_ug = mape(y_true, y_pred_ug)
    mae_ug = mae(y_true, y_pred_ug)

    print(Fore.YELLOW + "\n📊 Summary of Errors:")
    print(f"SSG → MAPE = {mape_ssg:.2f} %, MAE = {mae_ssg:.3f} kW/m²")
    print(f"UG  → MAPE = {mape_ug:.2f} %, MAE = {mae_ug:.3f} kW/m²")

    return y_pred_ssg, y_pred_ug, {
        "SSG": {"MAPE (%)": mape_ssg, "MAE": mae_ssg},
        "UG":  {"MAPE (%)": mape_ug, "MAE": mae_ug}
    }