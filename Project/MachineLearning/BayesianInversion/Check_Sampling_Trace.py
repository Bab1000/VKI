import numpy as np
import arviz as az # type: ignore
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde



# === Getting the samples ===
trace_path = f"/home/jpe/VKI/Project/MachineLearning/BayesianInversion/New_Simulations/Mass_Flow_based_testing/UniformG_MultiCond_Pstat={1500/100}/MassFlow={10}/Posteriors/Thinned_trace.nc"  
idata = az.from_netcdf(trace_path)

def plot_trace_and_pdf(idata, var_names):
    """
    Plot trace and KDE of posterior samples from ArviZ idata for given variables.

    Parameters:
    - idata: ArviZ InferenceData object
    - var_names: list of variable names to plot (e.g., ["GN", "GO"])
    """
    posterior = idata.posterior
    n_chains = posterior.sizes["chain"]
    colors = plt.cm.tab10.colors

    for name in var_names:
        if name not in posterior.data_vars:
            print(f"[WARN] Variable {name} not found in posterior. Skipping.")
            continue

        samples = posterior[name]  # (chain, draw)

        fig, axs = plt.subplots(1, 2, figsize=(12, 4))

        for c in range(n_chains):
            chain_samples = samples.sel(chain=c).values.flatten()

            # Trace plot
            axs[0].plot(chain_samples, label=f"Chain {c}", color=colors[c % len(colors)], alpha=0.7)

            # KDE
            kde = gaussian_kde(chain_samples)
            x_vals = np.linspace(np.min(chain_samples), np.max(chain_samples), 300)
            y_vals = kde(x_vals)
            axs[1].plot(x_vals, y_vals, color=colors[c % len(colors)])
            axs[1].fill_between(x_vals, y_vals, alpha=0.2, color=colors[c % len(colors)])

        axs[0].set_title(f"Trace - {name}")
        axs[0].set_ylabel("Value")

        axs[1].set_title(f"PDF - {name}")
        axs[1].set_ylabel("Density")

        plt.tight_layout()
        plt.show()

plot_trace_and_pdf(idata, ["GN", "GO"])

