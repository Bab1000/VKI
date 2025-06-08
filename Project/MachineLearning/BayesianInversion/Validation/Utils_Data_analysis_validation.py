from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns # type: ignore
from colorama import Fore
import arviz as az # type: ignore

def Mean_plots(sm_model,T_test,HF_test,HF_error,Pdyn_test,Pstat_test,trace_path,output_path):
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
        print(Fore.WHITE + f"---> [INFO] Computing mean HF : Pstat = {Pstat_test[i]/100} | Pdyn = {Pdyn_test[i]:.1f} | T = {T_test[i]:.1f} ...")

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

def MAP_plots(sm_model,samples,Pdyn_test,T_test,T_exp_list,HF_exp_list,HF_exp_err_list,output_path,std_pred=None,method="global",n_local_samples=5000,annotate=True):
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
        label="Experimental (± offset)"
    )

    for i in range(len(T_exp_list)):
        try:

            print(Fore.WHITE + f"---> [INFO] Computing MAP HF : Pdyn = {Pdyn_test[i]:.1f} | T = {T_test[i]:.1f} ...")

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
            map_hf = get_map(predictions,T_exp_list[i],HF_exp_list[i])
            min_pred = predictions.min()
            max_pred = predictions.max()
            mid_pred = (min_pred + max_pred) / 2
            err_pred = [[mid_pred - min_pred], [max_pred - mid_pred]]  # Asymmetric error bars

            plt.errorbar(
                T_test[i],
                map_hf,
                yerr=[[map_hf - min_pred], [max_pred - map_hf]],
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
    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

def get_map(predictions,T,HF_exp_val):
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

    output_path=f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case/Plots/HF_dist_T={T:.0f}.jpeg"

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
    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

    return map_val

def plot_joint_GN_GO(samples):
    """
    Plots the joint distribution of GN and GO from posterior samples.

    Parameters:
        samples: dict-like or ArviZ InferenceData object with "GN" and "GO"
        kind: "kde" or "scatter"
        bins: bin size (if kind='hexbin')
        show: whether to display the plot
        save_path: path to save the figure (if not None)
    """

    print(Fore.WHITE + f"---> [INFO] Computing GN-GO joint distribution ...")

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

    output_path=f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/Validation/Numerical_test_case/Plots/GN_GO_JointDist.jpeg"

    plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

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
        plt.savefig(global_output_path + f"/GN_GO_{targ_p/100}.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
    else:
        print("[WARNING] One or both of GN and GO arrays are empty.")
