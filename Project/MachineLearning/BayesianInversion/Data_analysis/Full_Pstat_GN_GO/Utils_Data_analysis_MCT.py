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

def MAP_plots(sm_model, samples, targ_p, MF_test, Pdyn_test, T_test,
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

    # === Setup for MAP points ===
    color_dict = {10: 'blue', 16: 'green', 20: 'red'}
    ecolor_dict = {10: 'lightblue', 16: 'lightgreen', 20: 'lightcoral'}
    marker_dict = {10: 's', 16: 'D', 20: '^'}  # square, diamond, triangle

    labels_done = set()

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

            mf_value = round(MF_test[i])
            for ref_mf in [10, 16, 20]:
                if abs(mf_value - ref_mf) <= 3:
                    label = f"MAP ± range (MF = {ref_mf} g/s)"
                    if label in labels_done:
                        label = None
                    else:
                        labels_done.add(label)

                    plt.errorbar(
                        T_test[i],
                        map_hf,
                        yerr=std_pred,
                        fmt=marker_dict[ref_mf],
                        color=color_dict[ref_mf],
                        ecolor=ecolor_dict[ref_mf],
                        elinewidth=2,
                        capsize=4,
                        markersize=7,
                        alpha=0.9,
                        label=label
                    )

                    #Optional annotation
                    #plt.annotate(f"{mf_value}g/s", (T_test[i], map_hf),
                    #             textcoords="offset points", xytext=(0, 5),
                    #             ha='center', fontsize=8)

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

    output_path_HF = global_output_path + f"/HF_distributions/FullCond_MAP_Pstat={targ_p/100}_T={T:.2f}.jpeg"

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

def plot_joint_GN_GO(samples, targ_p,global_output_path):
    """
    Plots the joint distribution of GN and GO from posterior samples.

    Parameters:
        samples: dict-like or ArviZ InferenceData object with "GN" and "GO"
        kind: "kde" or "scatter"
        bins: bin size (if kind='hexbin')
        show: whether to display the plot
        save_path: path to save the figure (if not None)
    """
    print(Fore.BLUE + f"[STEP] Computing joint probabilities of GN and GO for Pstats = {targ_p}/100 mbar")

    print(Fore.WHITE + f"---> [INFO] Computing GN-GO joint distribution : Pstat = {targ_p/100} ...")

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

    output_path=global_output_path + f"/GN_GO_Jointdistributions/GN_GO_JointDist={targ_p/100}.jpeg"

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

def plot_correlation(trace, Pstat, common_vars=('GN', 'GO'), kind='scatter',
                                   show_plots=True, save_dir=None):
    """
    Plots a correlation matrix for each test case in a PyMC trace.

    Parameters:
    - trace: arviz.InferenceData
        The trace from pm.sample() in InferenceData format.
    - common_vars: tuple of str
        Variable names that are shared across all test cases.
    - method: str
        Correlation method: 'pearson' or 'spearman'.
    - show_plots: bool
        Whether to display the plots with plt.show().
    - save_dir: str or None
        If provided, saves plots as PNG files in this directory.
    """

    print(Fore.BLUE + f"[STEP] Computing posteriors correlations for Pstat = {Pstat/100} mbar")

    all_vars = trace.posterior.data_vars.keys()

    # Find unique test case suffixes like "TC=1", "TC=2", etc.
    tc_suffixes = sorted(set(
        re.search(r'TC=\d+', var).group()
        for var in all_vars if 'TC=' in var
    ))

    # Filter common variables to those that actually exist
    common_present = [v for v in common_vars if v in all_vars]

    for tc in tc_suffixes:

        print(Fore.WHITE + f"---> [INFO] Test case {tc} analysis ...")

        # Variables specific to this test case
        tc_vars = [var for var in all_vars if var.endswith(f"_{tc}")]
        selected_vars = tc_vars + common_present

        if len(selected_vars) < 2:
            continue  # Skip if not enough variables to correlate

        az.plot_pair(
            trace,
            var_names=selected_vars,
            kind=kind,                   # 'scatter' or 'kde'
            marginals=True,
            divergences=False,
            backend_kwargs={"figsize": (10, 10)}
        )
        #plt.tight_layout()

        if save_dir:
            fname = f"{save_dir}/correlation_matrix_{tc}_Pstat={Pstat}.jpeg"
            plt.savefig(fname)

        if show_plots:
            plt.show()
        else:
            plt.close()

def sensitivity_gamma(sm_model, gamma_range=np.linspace(-4, 0, 25),
                      Pdyn_list=None, MF_list=None, T_list=None, target_p=None,
                      output_path=None):
    """
    Plot HF vs T with min–max error bars from gamma variation.
    Marker shape depends on mass flow (MF_list).
    Color: blue for gammaO varied, red for gammaN varied.

    Parameters:
    - sm_model: surrogate model with .predict_values(X)
    - gamma_range: range of gamma values to test
    - Pdyn_list: list of dynamic pressures
    - MF_list: list of mass flow rates [g/s] (used for markers)
    - T_list: list of temperatures
    - target_p: fixed static pressure (in Pa)
    - output_path: optional path to save figure
    """

    print(Fore.BLUE + f"[STEP] Gammas sensitivity analysis study")

    if any(v is None for v in [Pdyn_list, MF_list, T_list]):
        raise ValueError("Pdyn_list, MF_list, and T_list must all be provided")
    assert len(Pdyn_list) == len(MF_list) == len(T_list), "All input lists must have same length"

    marker_dict = {10: 'o', 16: 's', 20: '^'}
    color_varied = {
        'GO': 'red',  # GO fixed ⇒ GN varying ⇒ red
        'GN': 'blue'  # GN fixed ⇒ GO varying ⇒ blue
    }

    plt.figure(figsize=(8, 5))

    for fix_gamma in ['GO', 'GN']:
        print(Fore.WHITE + f"---> [INFO] Computing for fixed {fix_gamma} ...")
        for i in range(len(T_list)):
            T = T_list[i]
            Pdyn = Pdyn_list[i]
            MF = round(MF_list[i])
            hf_samples = []

            for gamma in gamma_range:
                GN, GO = (-4, gamma) if fix_gamma == 'GN' else (gamma, -4)
                X = np.array([[Pdyn, target_p, T, GN, GO]])
                HF = sm_model.predict_values(X)[0, 0] / 1000
                hf_samples.append(HF)

            hf_samples = np.array(hf_samples)
            hf_min, hf_max = np.min(hf_samples), np.max(hf_samples)
            mean_hf = np.mean(hf_samples)
            err_low = mean_hf - hf_min
            err_high = hf_max - mean_hf

            #marker = marker_dict.get(MF, 'x')
            color = color_varied[fix_gamma]
            label = None
            if i == 0:
                # label reflects the gamma that is VARYING
                label = f"G{'N' if fix_gamma == 'GO' else 'O'} varied"

            plt.errorbar(
                T, mean_hf,
                yerr=[[err_low], [err_high]],
                fmt="none",       #marker
                color=color,
                ecolor=color,
                alpha=0.4,
                capsize=4,
                markersize=6,
                label=label
            )

    plt.xlabel("Temperature [K]")
    plt.ylabel("Predicted Heat Flux [kW/m²]")
    plt.title(f"Gamma sensitivity analysis | Pstat = {target_p / 100:.0f} mbar")
    plt.grid(True, linestyle=':')

    legend_elements = [
        # Gamma variations (color only)
        mlines.Line2D([], [], color='red', marker='o', linestyle='None', label='GN varied | GO = -4'),
        mlines.Line2D([], [], color='blue', marker='o', linestyle='None', label='GO varied | GN = -4'),
        
        # Mass flow markers (shape only, black color)
        mlines.Line2D([], [], color='black', marker='o', linestyle='None', label='MF = 10 g/s'),
        mlines.Line2D([], [], color='black', marker='s', linestyle='None', label='MF = 16 g/s'),
        mlines.Line2D([], [], color='black', marker='^', linestyle='None', label='MF = 20 g/s'),
    ]

    plt.legend(handles=legend_elements, loc='upper left')

    plt.tight_layout()

    if output_path:
        fname = f"{output_path}/Gamma_sensitivity/Gamma_sensitivity_Pstat={target_p}.jpeg"
        plt.savefig(fname, format='jpeg', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

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

        plt.title("KDE of GN and GO")
        plt.xlabel("log10(Gamma)")
        plt.ylabel("Probability Density")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(global_output_path + f"/GN_GO/{targ_p/100}.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
    else:
        print("[WARNING] One or both of GN and GO arrays are empty.")

def plot_T_analysis(trace, T0_list, percent_error, global_output_path,targ_p, n_cases=5):
    """
    Plot posterior distributions of T and Gaussian priors computed from % uncertainty.

    Parameters:
    - trace: posterior samples (from pm.sample())
    - T0_list: list of nominal T values (≥ n_cases)
    - percent_error: relative error in %, used to compute σ such that 95%CI = ± ε%
    - n_cases: number of test cases to plot
    """

    colors = plt.cm.tab10(np.linspace(0, 1, n_cases))
    plt.figure(figsize=(10, 6))

    T0_array = np.array(T0_list)  # au cas où c'était une liste
    indices = np.random.choice(len(T0_array), size=n_cases, replace=False)

    it = 0
    for i in indices:
        key = f"T_TC={i}"
        if key not in trace:
            print(f"[WARNING] {key} not in trace.")
            continue

        T0 = T0_list[i]
        sigma = (percent_error) * T0 / 1.96  # convert % to std

        # Posterior
        post_samples = trace[key].values.flatten()
        x_post = np.linspace(post_samples.min(), post_samples.max(), 500)
        y_post = gaussian_kde(post_samples)(x_post)

        # Prior
        x_prior = np.linspace(T0 - 4*sigma, T0 + 4*sigma, 500)
        y_prior = norm.pdf(x_prior, loc=T0, scale=sigma)

        color = colors[it]
        label = f"Test case {i}"

        # Plot posterior
        plt.plot(x_post, y_post, color=color, label=f"{label} posterior")
        plt.fill_between(x_post, y_post, color=color, alpha=0.3)

        # Plot prior
        plt.plot(x_prior, y_prior, linestyle='--', color=color, alpha=0.6, label=f"{label} prior")

        it = it + 1

    plt.xlabel("T [K]")
    plt.ylabel("Probability Density")
    plt.title(f"Posterior vs Prior of T")
    plt.grid(True)
    legend_elements = [
        Line2D([], [], color='black', linestyle='-', label='Posterior'),
        Line2D([], [], color='black', linestyle='--', label='Prior'),
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.tight_layout()
    plt.savefig(global_output_path + f"/T_analysis/{targ_p/100}.jpeg", format='jpeg', dpi=300, bbox_inches='tight')

def plot_relative_errors_by_massflow(sm_model, samples, targ_p, T_test, HF_exp_list, MF_test, output_path):
    """
    Plots the relative error (%) between predicted and experimental heat flux
    using line plots grouped and colored by mass flow.
    Also computes and displays the weighted average relative error.
    """

    import matplotlib.pyplot as plt
    import numpy as np
    import os

    print(Fore.BLUE + "[STEP] Plotting relative errors (MAP vs Experimental) by mass flow")

    # === Setup mass flow (catégorique)
    mf_color_dict = {10: 'blue', 16: 'green', 20: 'red'}
    mf_marker_dict = {10: 's', 16: 'D', 20: '^'}

    grouped_data = {10: [], 16: [], 20: []}

    total_true_hf = 0
    total_abs_error = 0

    for i in range(len(T_test)):
        try:
            print(Fore.RED + f"---> [INFO] Colputing for TC={i} ...")
            # === Extraction des échantillons
            Pdyn_samples = samples[f"Pdyn_TC={i}"].values.flatten()
            Pstat_samples = samples[f"Pstat_TC={i}"].values.flatten()
            T_samples = samples[f"T_TC={i}"].values.flatten()
            GN_samples = samples["GN"].values.flatten()
            GO_samples = samples["GO"].values.flatten()

            # === Calcul MAP
            min_len = min(len(Pdyn_samples), len(Pstat_samples), len(T_samples), len(GN_samples), len(GO_samples))
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
            hf_exp = HF_exp_list[i]
            map_pred = np.median(predictions)

            rel_error = abs(map_pred - hf_exp) / hf_exp
            total_abs_error += abs(map_pred - hf_exp) / hf_exp

            total_true_hf += hf_exp

            # Identifier mass flow (arrondi proche)
            mf_value = round(MF_test[i])
            closest_mf = min(grouped_data.keys(), key=lambda x: abs(x - mf_value))

            grouped_data[closest_mf].append((T_test[i], rel_error))

        except Exception as e:
            print(Fore.RED + f"[ERROR] Failed to compute relative error for point {i}: {e}")
            continue

    # === Plot
    plt.figure(figsize=(8, 5))
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    for mf, values in grouped_data.items():
        if values:
            values.sort()
            temps, errors = zip(*values)
            temps = np.array(temps)
            errors = np.abs(np.array(errors))

            plt.plot(
                temps,
                errors,
                label=f"MF = {mf} g/s",
                color=mf_color_dict[mf],
                marker=mf_marker_dict[mf],
                linewidth=2,
                markersize=6
            )

    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.xlabel("T [K]")
    plt.ylabel("Relative Error (%)")
    plt.title(f"Relative Error (MAP vs Experimental) @ Pstat = {int(targ_p)} Pa")

    # === Affichage de l’erreur moyenne pondérée
    if total_true_hf > 0:
        avg_weighted_error = 100 * total_abs_error / len(T_test)
        print(Fore.MAGENTA + f"[MODEL EVAL] Estimated avg. relative error: {avg_weighted_error:.2f}%")

        plt.text(
            0.98, 0.02,
            f"Weighted avg. error = {avg_weighted_error:.2f}%",
            transform=plt.gca().transAxes,
            fontsize=10,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray')
        )
    else:
        print(Fore.YELLOW + "[WARNING] Cannot compute average model error (division by 0)")

    plt.legend()
    plt.tight_layout()

    # === Sauvegarde
    os.makedirs(os.path.join(output_path, "Errors"), exist_ok=True)
    save_path = os.path.join(output_path, f"Errors/relative_error_P{int(targ_p)}Pa.jpeg")
    plt.savefig(save_path, dpi=300)
    print(Fore.GREEN + f"[SUCCESS] Relative error plot saved to: {save_path}")
