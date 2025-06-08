from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns # type: ignore
from colorama import Fore
import arviz as az # type: ignore
import sys
import pandas as pd
from smt.surrogate_models import KRG


def Mean_plots(sm_model,targ_p,targ_mf,T_test,HF_test,HF_error,Pdyn_test,Pstat_test,trace_path,output_path):
    """
    Plot experimental data and surrogate model predictions using posterior samples.

    Parameters:
        sm_model: the surrogate model (e.g., Kriging)
        T_test, HF_test, HF_error: arrays of experimental data (converted to kW/m²)
        Pdyn_test, Pstat_test: arrays of inputs
        trace_path: path to .nc file with posterior trace
        output_path: where to save the figure (jpeg)
    """

    # --- Load trace
    idata = az.from_netcdf(trace_path)
    samples = idata.posterior.stack(sample=("chain", "draw"))

    plt.figure(figsize=(8, 6))

    # Plot experimental data with uncertainty
    plt.errorbar(
        T_test, HF_test,
        yerr=HF_error,
        fmt='o',
        color='black',
        ecolor='gray',
        elinewidth=1.5,
        capsize=4,
        label="Experimental data (± offset)"
    )

    for i in range(len(Pdyn_test)):
        print(Fore.WHITE + f"---> [INFO] Computing mean HF : Pstat = {targ_p/100} | MF = {targ_mf:.1f} | Pdyn = {Pdyn_test[i]:.1f} | T = {T_test[i]:.1f} ...")

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

        # Sample posterior
        n_points = 10000
        min_len = min(len(Pdyn_samples), len(Pstat_samples), len(T_samples), len(GN_samples), len(GO_samples))
        if min_len < n_points:
            raise ValueError(f"Not enough samples: requested {n_points}, but got {min_len}")

        indices = np.random.choice(min_len, size=n_points, replace=False)

        predictions = []
        input_pairs = []

        for idx in indices:
            XV = np.array([
                Pdyn_samples[idx],
                Pstat_samples[idx],
                T_samples[idx],
                GN_samples[idx],
                GO_samples[idx]
            ]).reshape(1, -1)

            pred = sm_model.predict_values(XV)[0, 0] / 1000
            predictions.append(pred)
            input_pairs.append((GN_samples[idx], GO_samples[idx]))

        predictions = np.array(predictions)
        mean_pred = predictions.mean()
        std_pred = predictions.std() * 1.96  # 95% confidence interval

        idx_closest = np.abs(predictions - mean_pred).argmin()
        gn_mean = 10 ** input_pairs[idx_closest][0]
        go_mean = 10 ** input_pairs[idx_closest][1]

        plt.errorbar(
            T_test[i],
            mean_pred,
            yerr=std_pred,
            fmt='o',
            color='tab:blue',
            ecolor='lightblue',
            elinewidth=2,
            capsize=4,
            label="Predicted (mean ± 95% CI)" if i == 0 else ""
        )

        plt.text(
            T_test[i] + 50, mean_pred + 200,
            f"GN={gn_mean:.5f}\nGO={go_mean:.5f}",
            fontsize=8,
            color="blue",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="blue", facecolor="white", alpha=0.7)
        )

    plt.legend()
    plt.title("Experimental vs Mean Predicted Heat Flux")
    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

def MAP_plots(sm_model,samples,targ_p,targ_mf,Pdyn_test,T_test,T_exp_list,HF_exp_list,HF_exp_err_list,output_path,global_output_path,std_pred=None,method="global",n_local_samples=5000,annotate=True):
    """
    Loops over all posterior samples and plots MAP predictions vs experimental data.

    Parameters:
        sm_model: surrogate model
        samples: ArviZ stacked posterior samples
        T_exp_list, HF_exp_list, HF_exp_err_list: lists of experimental values (same length)
        output_path: where to save the global MAP plot
        std_pred: global std dev (only needed for method='global')
        method: 'global' or 'local'
        n_local_samples: used if method='local'
        annotate: show GN/GO on plot
    """

    plt.figure(figsize=(8, 6))

    # Plot all experimental points
    plt.errorbar(
        T_exp_list, HF_exp_list,
        yerr=HF_exp_err_list,
        fmt='o',
        color='black',
        ecolor='gray',
        elinewidth=1.5,
        capsize=4,
        label="Experimental (± error)"
    )

    for i in range(len(T_exp_list)):
        try:

            print(Fore.WHITE + f"---> [INFO] Computing MAP HF : Pstat = {targ_p/100} | MF = {targ_mf:.1f} | Pdyn = {Pdyn_test[i]:.1f} | T = {T_test[i]:.1f} ...")

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

            # Sample posterior
            n_points = 10000
            min_len = min(len(Pdyn_samples), len(Pstat_samples), len(T_samples), len(GN_samples), len(GO_samples))
            if min_len < n_points:
                raise ValueError(f"Not enough samples: requested {n_points}, but got {min_len}")

            predictions = []
            input_pairs = []

            for idx in range(len(Pdyn_samples)):
                XV = np.array([
                    Pdyn_samples[idx],
                    Pstat_samples[idx],
                    T_samples[idx],
                    GN_samples[idx],
                    GO_samples[idx]
                ]).reshape(1, -1)

                pred = sm_model.predict_values(XV)[0, 0] / 1000
                predictions.append(pred)
                input_pairs.append((GN_samples[idx], GO_samples[idx]))

            predictions = np.array(predictions)
            map_hf = get_map(predictions,HF_exp_list[i],global_output_path,T_test[i],targ_p,targ_mf)
            std_pred = predictions.std() * 1.96

            plt.errorbar(
                T_test[i],
                map_hf,
                yerr=std_pred,
                fmt='s',
                color='red',
                ecolor='lightblue',
                elinewidth=2,
                capsize=4,
                label="MAP ± range" if i == 0 else ""
            )

        except Exception as e:
            print(Fore.RED + f"[ERROR] MAP computation failed for point {i}: {e}")

    plt.legend()
    plt.title("Experimental vs MAP Heat Flux")
    plt.xlabel("T [K]")
    plt.ylabel("HF [kW/m2]")
    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

def get_map(predictions,HF_exp_val,global_output_path,T,targ_p,targ_mf):
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

    output_path_HF = global_output_path + f"/HF_distributions/FullCond_MAP_Pstat={targ_p/100}_mf={targ_mf:.0f}_T={T:.2f}.jpeg"

    # Plot
    plt.figure(figsize=(8, 4))
    plt.plot(x_vals, y_vals, color='tab:blue')
    plt.axvline(map_val, color='red', linestyle='--', linewidth=2, label=f"MAP = {map_val:.2f} kW/m²")
    plt.axvline(HF_exp_val, color='green', linestyle='--', linewidth=2, label=f"Exp. value = {HF_exp_val:.2f} kW/m²")
    plt.title("Posterior HF Distribution")
    plt.xlabel("Predicted Heat Flux (kW/m²)")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path_HF, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

    return map_val

def plot_joint_GN_GO(samples, targ_p, targ_mf,global_output_path):
    """
    Plots the joint distribution of GN and GO from posterior samples.

    Parameters:
        samples: dict-like or ArviZ InferenceData object with "GN" and "GO"
        kind: "kde" or "scatter"
        bins: bin size (if kind='hexbin')
        show: whether to display the plot
        save_path: path to save the figure (if not None)
    """

    print(Fore.WHITE + f"---> [INFO] Computing GN-GO joint distribution : Pstat = {targ_p/100} | MF = {targ_mf:.1f} ...")

    try:
        GN = samples["GN"].values.flatten()
        GO = samples["GO"].values.flatten()
    except Exception as e:
        print(f"[ERROR] Cannot extract GN/GO from samples: {e}")
        return

    # Create 2D KDE
    values = np.vstack([GN, GO])
    kde = gaussian_kde(values)

    # Grid setup
    x = np.linspace(-4, 0, 200)
    y = np.linspace(-4, 0, 200)
    X, Y = np.meshgrid(x, y)
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = kde(positions).reshape(X.shape)

    # Plot
    plt.figure()
    contour = plt.contourf(X, Y, Z, cmap="viridis")
    plt.colorbar(contour, label="Density")

    plt.xlim(-4, 0)
    plt.ylim(-4, 0)
    plt.xlabel("GN (log10)")
    plt.ylabel("GO (log10)")
    plt.title("KDE Contour of GN and GO (gaussian_kde)")
    plt.tight_layout()

    output_path=global_output_path + f"/GN_GO_distributions/GN_GO_JointDist={targ_p/100}_mf={targ_mf:.0f}.jpeg"

    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

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
