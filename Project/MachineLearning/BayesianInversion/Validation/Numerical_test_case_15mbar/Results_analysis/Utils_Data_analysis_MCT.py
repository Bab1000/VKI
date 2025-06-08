from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde,norm
import seaborn as sns # type: ignore
from colorama import Fore
import arviz as az # type: ignore
import sys
import pandas as pd
from smt.surrogate_models import KRG
import re
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
import os


def MAP_plots(sm_model, samples, targ_p, Pdyn_test, T_test,
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
        elinewidth=1.8,
        capsize=4,
        markersize=7,
        label="Experimental (± error)"
    )

    for i in range(len(T_exp_list)):
        try:
            print(Fore.WHITE + f"---> [INFO] Computing MAP HF : Pstat = {targ_p/100} | Pdyn = {Pdyn_test[i]:.1f} | T = {T_test[i]:.1f} ...")

            try:
                Pdyn_samples = samples[f"Pdyn_TC={i}"].values.flatten()
                Pstat_samples = samples[f"Pstat_TC={i}"].values.flatten()
                T_samples = samples[f"T_TC={i}"].values.flatten()
                GN_samples = samples["GN"].values.flatten()
                GO_samples = samples["GO"].values.flatten()
            except KeyError as e:
                print(Fore.RED + f"[ERROR] Missing variable in trace: {e}")
                print(Fore.YELLOW + f"Available variables: {list(samples.keys())}")
                continue

            n_points = 10000
            min_len = min(len(Pdyn_samples), len(Pstat_samples), len(T_samples), len(GN_samples), len(GO_samples))
            if min_len < n_points:
                raise ValueError(f"Not enough samples: requested {n_points}, but got {min_len}")

            predictions = []
            for idx in range(min_len):
                XV = np.array([
                    Pdyn_samples[idx],
                    Pstat_samples[idx],
                    T_samples[idx],
                    GN_samples[idx],
                    GO_samples[idx]
                ]).reshape(1, -1)

                pred = sm_model.predict_values(XV)[0, 0] / 1000
                predictions.append(pred)

            predictions = np.array(predictions)
            map_hf = get_map(predictions, HF_exp_list[i], global_output_path, T_test[i], targ_p)
            std_pred = predictions.std() * 1.96

            plt.errorbar(
                T_test[i],
                map_hf,
                yerr=std_pred,
                elinewidth=2.2,
                capsize=5,
                markersize=8,
                alpha=0.9,
                color="#4A90E2",
                marker='s',
                label='Predicted HF' if i == 0 else ""
            )

        except Exception as e:
            print(Fore.RED + f"[ERROR] MAP computation failed for point {i}: {e}")

    # === Custom axes formatting ===
    plt.xlabel("T [K]", fontsize=16)
    plt.ylabel("HF [kW/m²]", fontsize=16)
    plt.tick_params(axis='both', labelsize=14, width=2, length=6, direction='inout')
    plt.legend(fontsize=13)
    plt.tight_layout()
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

    output_path_HF = global_output_path + f"/HF_distributions/FullCond_MAP_Pstat={targ_p/100}_T={T:.2f}.jpeg"

    # Plot
    plt.figure(figsize=(8, 4))
    plt.plot(x_vals, y_vals, color='tab:blue')
    plt.axvline(map_val, color='red', linestyle='--', linewidth=2, label=f"MAP = {map_val:.2f} kW/m²")
    plt.axvline(HF_exp_val, color='green', linestyle='--', linewidth=2, label=f"Exp. value = {HF_exp_val:.2f} kW/m²")
    plt.xlabel("Predicted Heat Flux (kW/m²)")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(output_path_HF, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

    return map_val

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

def GN_GO_plotting(samples,targ_p,global_output_path):
    try:
        GN = samples["GN"].values.flatten()
        GO = samples["GO"].values.flatten()
    except Exception as e:
        print(f"[ERROR] Cannot extract GN/GO from samples: {e}")

    if GN.size > 0 and GO.size > 0:
        # Create 1D KDEs
        kde_gn = gaussian_kde(GN)
        kde_go = gaussian_kde(GO)

        # Create evaluation range
        min_val = min(GN.min(), GO.min())
        max_val = max(GN.max(), GO.max())
        x_grid = np.linspace(min_val, max_val, 1000)

        # Evaluate KDEs
        y_gn = kde_gn(x_grid)
        y_go = kde_go(x_grid)

        # Define uniform prior over [-4, 0]
        min_gammaN = -4
        delta_gammaN = 4

        # Support of the prior
        x_prior = [min_gammaN, min_gammaN + delta_gammaN]  # [-4, 0]

        # Constant density value for uniform distribution
        y_prior = [1 / delta_gammaN, 1 / delta_gammaN]  # [0.25, 0.25]

        # Plot both KDEs
        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(x_grid, y_gn, color="#4A90E2", label="GN")  # bleu pastel
        plt.fill_between(x_grid, y_gn, color="#4A90E2", alpha=0.3)

        plt.plot(x_grid, y_go, color="#7ED6A5", label="GO")  # vert pastel
        plt.fill_between(x_grid, y_go, color="#7ED6A5", alpha=0.3)

        plt.plot(x_prior, y_prior, color='red', linewidth=1.5, linestyle='--', label='Prior')
        plt.xlabel("log10(Gamma)")
        plt.ylabel("Probability Density")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(global_output_path + f"/GN_GO/{targ_p/100}.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
    else:
        print("[WARNING] One or both of GN and GO arrays are empty.")

def GN_GO_separate_plotting(samples, targ_p, global_output_path):
    try:
        GN = samples["GN"].values.flatten()
        GO = samples["GO"].values.flatten()
    except Exception as e:
        print(f"[ERROR] Cannot extract GN/GO from samples: {e}")
        return

    if GN.size > 0 and GO.size > 0:
        from scipy.stats import gaussian_kde
        import numpy as np
        import matplotlib.pyplot as plt
        import os

        # Create KDEs
        kde_gn = gaussian_kde(GN)
        kde_go = gaussian_kde(GO)

        # Shared evaluation range
        min_val = min(GN.min(), GO.min())
        max_val = max(GN.max(), GO.max())
        x_grid = np.linspace(min_val, max_val, 1000)

        # Evaluate KDEs
        y_gn = kde_gn(x_grid)
        y_go = kde_go(x_grid)

        # Define prior (uniform on [-4, 0])
        min_gammaN = -4
        delta_gammaN = 4
        x_prior = [min_gammaN, min_gammaN + delta_gammaN]
        y_prior = [1 / delta_gammaN, 1 / delta_gammaN]  # [0.25, 0.25]

        # Create save directory
        save_dir = os.path.join(global_output_path, "GN_GO")
        os.makedirs(save_dir, exist_ok=True)

        # === GN Figure ===
        plt.figure(figsize=(8, 5))
        plt.plot(x_grid, y_gn, color="#4A90E2", label="GN", linewidth=2.5)
        plt.fill_between(x_grid, y_gn, color="#4A90E2", alpha=0.3)
        plt.plot(x_prior, y_prior, color='red', linewidth=2, linestyle='--', label='Prior')

        plt.xlabel("log10(Gamma)", fontsize=14)
        plt.ylabel("Probability Density", fontsize=14)

        # ➕ Augmenter la taille des ticks
        plt.tick_params(axis='both', labelsize=12, width=2, length=6, direction='inout')

        plt.grid(True)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"GN_validation.jpeg"), format='jpeg', dpi=300, bbox_inches='tight')
        plt.close()


        # === GO Figure ===
        plt.figure(figsize=(8, 5))
        plt.plot(x_grid, y_go, color="#7ED6A5", label="GO", linewidth=2.5)
        plt.fill_between(x_grid, y_go, color="#7ED6A5", alpha=0.3)
        plt.plot(x_prior, y_prior, color='red', linewidth=2, linestyle='--', label='Prior')

        plt.xlabel("log10(Gamma)", fontsize=14)
        plt.ylabel("Probability Density", fontsize=14)

        # ➕ Augmenter la taille des ticks
        plt.tick_params(axis='both', labelsize=12, width=2, length=6, direction='inout')

        plt.grid(True)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"GO_validation.jpeg"), format='jpeg', dpi=300, bbox_inches='tight')
        plt.close()


    else:
        print("[WARNING] One or both of GN and GO arrays are empty.")
