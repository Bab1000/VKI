import pytensor  # type: ignore
pytensor.config.mode = "FAST_COMPILE"
import pymc as pm # type: ignore
import numpy as np
import arviz as az # type: ignore
import matplotlib.pyplot as plt
from colorama import Fore, Style, init
import sys
import warnings
from collections import defaultdict
import re
import scipy.stats as stats
import matplotlib.cm as cm
import os

# Initialize colorama
init(autoreset=True)

# ==============================================================================================
# ------------------------
# | Bayesian model class |
# ------------------------

class BayesianInversion:

    def __init__(self, y, forward_model):
        print(Fore.BLUE + "[STEP] Initializing the Bayesian model")
        try:
            
            print(Fore.WHITE + "---> [INFO] Building initialization ...")

            self.y = y                         # Observed features
            self.priors = []                   # Prior distributions
            self.observed = []                 # Observed distributions
            self.model = None                  # Bayesian model
            self.forward_model = forward_model # Forward surrogate model
            self.posteriors = None             # Posterior results

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Initialization failed: {e}")
            sys.exit(1)
            raise

    def build_model(self):
        print(Fore.WHITE + "---> [INFO] Building the Bayesian framework ...")
        try:
            with pm.Model() as model:

                # ===========================================================================================
                
                # Definition of prior distributions
                # ---------------------------------

                try:
                    for dist in self.priors:

                        # Gathering distribution infos
                        name = dist["name"]
                        dist_type = dist["type"]
                        mu = dist.get("mu")
                        uncertainty = dist.get("uncertainty")
                        lower = dist.get("lower")
                        upper = dist.get("upper")

                        # Construction of normal distribution
                        if dist_type == "normal" and (mu is not None and uncertainty is not None):
                            # Computation of standard deviation (95%)
                            sigma = max(uncertainty / 1.96, 1e-3)
                            # Building normal distribution
                            var = pm.Normal(name, mu=mu, sigma=sigma)

                        # Construction of normal distribution
                        elif dist_type == "uniform" and (lower is not None and upper is not None):
                            # Building uniform distribution
                            var = pm.Uniform(name, lower=lower, upper=upper)

                        else:
                            raise ValueError(f"Unsupported distribution: {dist_type}")

                except Exception as e:
                    print(Fore.RED + f"---> [ERROR] Prior {name} distribution was not added to the model: {e}")
                    sys.exit(1)

                # ===========================================================================================

                # Forward model definition (Stagline surrogate model)
                # --------------------------------------------------

                try:
                    # Check if this is a vectorized case (multiple outputs --> Check if 1D vector with more than one entry)
                    if len(np.shape(self.y)) == 1 and len(self.y) > 1:
                        print(Fore.WHITE + "---> [INFO] Detected vectorized setup (multiple experiment)")

                        # Building the instance linking pyMC with SMT
                        Stagline_SM = SMStaglineVectorized(self.forward_model)

                        # Check for grouped variables by their TC number or global variables
                        grouped_vars = defaultdict(list)
                        global_vars = defaultdict(list)

                        # Patter of test case number for multiple experiments recognition
                        pattern = r"TC=(\d+)"

                        for name in model.named_vars:
                            # Search for the pattern in the distribution's name
                            match = re.search(pattern, name)
                            if match:
                                # Gathering test case number
                                tc_number = int(match.group(1))
                                # Grouping distribution by TC number 
                                grouped_vars[tc_number].append(name)
                            else:
                                # Taking the distribution as global one if no TC number is recognized
                                global_vars["global"].append(name)


                        # Sort TC keys to keep consistent order
                        sorted_tc = sorted(grouped_vars.keys())
                        print(Fore.WHITE + f"---> [INFO] Grouped {len(sorted_tc)} Test Cases sets for vectorized modeling")
                        
                        # Building the input matrix for surrogate model
                        input_matrix = []
                        # Gather infos for each test case
                        for tc in sorted_tc:
                            row = []
                            for var_name in grouped_vars[tc]:
                                # Gather each distribution of the test case
                                var = model.named_vars.get(var_name)
                                if var is None:
                                    raise ValueError(f"Prior '{var_name}' not found in model.")
                                # Appending the values in the matrix
                                vec = pt.flatten(var)[0].dimshuffle(())
                                row.append(vec)


                            # Append global variables 
                            for var_name in global_vars["global"]:
                                # Gathering each global distribution 
                                var = model.named_vars.get(var_name)
                                # Appending the values in the matrix
                                vec = pt.flatten(var)[0].dimshuffle(())
                                row.append(vec)
                            
                            # Creating the global input matrix for the surrogate model
                            input_matrix.append(pt.stack(row))
                        
                        # Creating a pytensor input for the SMTStaglineVectorized class
                        SM_input = pt.stack(input_matrix, axis=0)  # shape: (n_exp, n_inputs)
                        # Prediction of the output of the model
                        SM_prediction = Stagline_SM(SM_input)

                    else:
                        # Non-vectorized case
                        print(Fore.WHITE + "---> [INFO] Detected scalar setup (single experiment)")
                        Stagline_SM = SMStagline(self.forward_model)

                        # Gathering the prior distributions
                        prior_var_list = []
                        for name in model.named_vars:
                            # Gathering the values for each prior distributions
                            var = model.named_vars.get(name)
                            if var is None:
                                raise ValueError(f"Prior distribution '{name}' not found in model.")
                            # Building the input vector for the SMT surrogate model
                            prior_var_list.append(pt.flatten(var))

                        # Creating a pytensor input for the SMTStaglineVectorized class
                        SM_input = pt.concatenate(prior_var_list, axis=0)  # shape: (n_inputs,)
                        # Prediction of the output of the model
                        SM_prediction = Stagline_SM(SM_input)


                except Exception as e:
                    print(Fore.RED + f"---> [ERROR] Failed to apply surrogate model: {e}")
                    sys.exit(1)


                # ===========================================================================================

                # Definition of observed distributions
                # ------------------------------------

                try:
                    counter = 0
                    if len(np.shape(self.y)) == 1 and len(self.y) > 1: 
                        for dist in self.observed:

                            # Gathering distribution infos
                            name = dist["name"]
                            dist_type = dist["type"]
                            mu = dist.get("mu")
                            uncertainty = dist.get("uncertainty")
                            lower = dist.get("lower")
                            upper = dist.get("upper")

                            # Construction of normal distribution
                            if dist_type == "normal" and (mu is not None and uncertainty is not None):
                                # Computation of standard deviation (95%)
                                sigma = max(uncertainty / 1.96, 1e-2)
                                # Uses Stagline SM for the observed quantities
                                var = pm.Normal(name, mu=SM_prediction[counter], sigma=sigma, observed = self.y[counter])

                            # Construction of normal distribution
                            elif dist_type == "uniform" and (lower is not None and upper is not None):
                                # Bulding uniform distribution
                                var = pm.Uniform(name, lower=lower, upper=upper, observed = self.y[counter])

                            else:
                                raise ValueError(f"Unsupported distribution: {dist_type}")
                    
                    else:
                        for dist in self.observed:

                            # Gathering distribution infos
                            name = dist["name"]
                            dist_type = dist["type"]
                            mu = dist.get("mu", 0)
                            uncertainty = dist.get("uncertainty", 1)

                            # Construction of normal distribution
                            if dist_type == "normal":
                                # Computation of standard deviation (95%)
                                sigma = max(uncertainty / 1.96, 1e-2)
                                # Uses Stagline SM for the observed quantities
                                var = pm.Normal(name, mu=SM_prediction, sigma=sigma, observed = self.y)

                            # Construction of normal distribution
                            elif dist_type == "uniform":
                                # Lower bound
                                lower = mu 
                                # Upper bound
                                upper = uncertainty 
                                # Building uniform distribution
                                var = pm.Uniform(name, lower=lower, upper=upper, observed = self.y)

                            else:
                                raise ValueError(f"Unsupported distribution: {dist_type}")


                except Exception as e:
                    print(Fore.RED + f"---> [ERROR] Observed {name} distribution was not added to the model: {e}")
                    sys.exit(1)

                # Construction of the Bayesian model
                self.model = model

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Model building failed: {e}")
            sys.exit(1)
    
    def add_priors(self, name, distribution="normal", mu=None, uncertainty=None, lower = None, upper = None):
        try:
            distribution = distribution.lower()

            if distribution == "normal" and (mu is None or uncertainty is None):
                raise ValueError("mu and uncertainty must be provided for a Normal distribution.")
            if distribution == "uniform" and (lower is None or upper is None):
                raise ValueError("uncertainty must be provided for a Uniform distribution.")

            # On stocke la définition de la distribution dans une liste
            self.priors.append({
                "name": name,
                "type": distribution,
                "mu": mu,
                "uncertainty": uncertainty,
                'lower': lower,
                "upper": upper
            })

            print(Fore.WHITE + f"---> [INFO] Prior {distribution.capitalize()} distribution registered : {name}")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register distribution '{name}': {e}")
            sys.exit(1)

    def add_observed(self, name, distribution="normal", mu=None, uncertainty=None, lower = None, upper = None):
        try:
            distribution = distribution.lower()

            if distribution == "normal" and (mu is None or uncertainty is None):
                raise ValueError("mu and uncertainty must be provided for a Normal distribution.")
            if distribution == "uniform" and uncertainty is None:
                raise ValueError("uncertainty must be provided for a Uniform distribution.")

            # On stocke la définition de la distribution dans une liste
            self.observed.append({
                "name": name,
                "type": distribution,
                "mu": mu,
                "uncertainty": uncertainty,
                'lower': lower,
                "upper": upper
            })

            print(Fore.WHITE + f"---> [INFO] Observed {distribution.capitalize()} distribution registered : {name}")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register distribution '{name}': {e}")
            sys.exit(1)

    def run_inference(self, draws = 2000, tune = 1000, chains = 4, cores = 4, restart_path = None, save_path = None):
        print(Fore.BLUE + "[STEP] Running Bayesian inference")
        try:
            if self.model is None:
                print(Fore.WHITE + "---> [INFO] Building the model instance ...")
                self.build_model()
                print(Fore.GREEN + "---> [SUCCESS] Model instance successfully constructed")

            # Starts the Inversion process using Markov Chain Monte Carlo process (MCMC)

            # Tune (Warm-up parameter):
            #       - MCMC has to "warm up" (adapt its internal settings like step size)
            #       - These early computations will not be included in the final results
            # Draws:
            #       - MCMC Sampling points

            warnings.filterwarnings("ignore", category=RuntimeWarning)
            np.seterr(over='ignore', invalid='ignore')
            
            # === Restart simulation from previous one ===
            start = None
            if restart_path:
                # Check for the restart file 
                path_restart_file = os.path.join(restart_path,"Restart_trace.nc")
                print(Fore.WHITE + f"---> [INFO] Resuming from saved trace: {path_restart_file}")
                previous = az.from_netcdf(path_restart_file)
                start = {
                    # Computing the mean values based on each draw of each chain
                    var: previous.posterior[var].mean(dim=["chain", "draw"]).values
                    for var in previous.posterior
                }

            # === Running MCMC ===

            with self.model:
                # === Creating sampling arguments ===
                step = pm.Metropolis()
                sampling_args = {
                    "draws": draws,                   # Number of time steps per parameters
                    "tune": tune,                     # Number of Warm-up or burn-in steps per chains
                    "chains": chains,                 # Number of chains
                    "cores": cores,                   # Number of cores
                    "return_inferencedata": True,
                    "step": step                      # Step to chose the next values
                    }
                
                if start is not None:
                    sampling_args["initvals"] = start    # Starting from a previous simulation option

                # Starting sampling
                print(Fore.WHITE + "---> [INFO] Starting Bayesian inversion (MCMC sampling) ...")
                # MCMC process
                self.posteriors = pm.sample(**sampling_args)

            print(Fore.GREEN + "---> [SUCCESS] Computation successfully completed !")

            # === Autosave for posteriors ===
            self.save_trace(save_path=save_path)

            # === Premilinary data analysis ===
            map_values, R_hat_all = self.priliminary_analysis()

            return self.posteriors, map_values, R_hat_all
        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Inference failed: {e}")

    def save_trace(self, save_path=None, max_points=10000, save_full=False, save_means=True):
        """
        Save the posterior trace in different formats:
        - Thinned version for plotting (reduced size)
        - Posterior means only (lightweight restart)
        - Optionally full trace (for full restart)

        Parameters:
        - save_path (str): base path to save (.nc, .csv will be appended)
        - max_points (int): max number of posterior samples to keep in thinned version
        - save_full (bool): whether to save the full trace or not
        - save_means (bool): whether to save posterior means as restart file
        """
        if save_path is None or self.posteriors is None:
            print(Fore.YELLOW + "---> [WARNING] No trace or path provided, skipping save.")
            return

        # Create directory if needed
        os.makedirs(save_path, exist_ok=True)
        idata = self.posteriors
        posterior = idata.posterior

        # === Random thinning for plotting purposes ===
        try:
            print(Fore.WHITE + "---> [INFO] Applying random thinning of posterior samples ...")
            # Gathering data
            n_chains = posterior.sizes["chain"]
            n_draws = posterior.sizes["draw"]
            total_samples = n_chains * n_draws
            # Maximum points to keep
            max_keep = int(min(max_points, 0.9*total_samples))
            print(Fore.WHITE + f"---> [INFO] Keeping {max_keep} points in total ...")
            # Number of points to keep in each chains
            draws_to_keep_per_chain = min(n_draws, max(max_keep // n_chains, 1))
            print(Fore.WHITE + f"---> [INFO] Keeping {draws_to_keep_per_chain} points in each of the {posterior.sizes['chain']} chains ...")

            # Final safety check
            if draws_to_keep_per_chain > n_draws:
                raise ValueError(f"Cannot keep {draws_to_keep_per_chain} draws per chain, only {n_draws} available.")
            
            # Starting the sampling of the chains
            thinned_data = {}
            for var in posterior.data_vars:
                var_data = posterior[var].values  # shape: (chain, draw)
                chain_list = []

                for chain_id in range(posterior.sizes["chain"]):
                    # Randomly chosing the points to keep
                    draw_indices = np.random.choice(n_draws, size=draws_to_keep_per_chain, replace=False)
                    draw_indices.sort()
                    # Gathering the thinned chain
                    thinned_chain = var_data[chain_id, draw_indices]  # shape: (n_draws,)
                    chain_list.append(thinned_chain)

                # Reconstructing the structure of the posteriors chains
                thinned_array = np.stack(chain_list, axis=0)  # shape: (chain, n_draws)
                thinned_data[var] = thinned_array

            # Build InferenceData to save it as .nc
            idata_thinned = az.from_dict(posterior=thinned_data)
            # Creating the full path of the file
            thinned_path = os.path.join(save_path, "Thinned_trace.nc")
            # Saving the file
            az.to_netcdf(idata_thinned, thinned_path)
            print(Fore.WHITE + f"---> [INFO] Thinned posterior saved to: {thinned_path}")

        except Exception as e:
            print(Fore.RED + f"[ERROR] Failed to save the thinned trace: {e}")

        # === Save posterior means only (for restart) ===
        if save_means:
            try:
                print(Fore.WHITE + "---> [INFO] Saving simulation trace for further restart purposes")
                # Computing the mean value of each prior
                mean_data = {}
                for var in posterior.data_vars:
                    # Get the mean value
                    mean_val = idata.posterior[var].values.mean()
                    # Store as a numpy array with shape (1, 1)
                    mean_data[var] = np.array([[mean_val]])  # shape: (chain=1, draw=1)

                # Creating inference data to be able to save it as .nc 
                idata_means = az.from_dict(posterior=mean_data)
                # Creating the full path of the file
                mean_path = os.path.join(save_path,"Restart_trace.nc")
                # Saving the file
                az.to_netcdf(idata_means, mean_path)
                print(Fore.WHITE + f"---> [INFO] Posterior means saved to: {mean_path}")
            except Exception as e:
                print(Fore.RED + f"[ERROR] Failed to save posterior means: {e}")

        # === Optionally save full trace ===
        if save_full:
            # Creating the full path of the file
            full_path = os.path.join(save_path,"Full_trace.nc")
            # Saving the file
            az.to_netcdf(idata, full_path)
            print(Fore.WHITE + f"---> [INFO] Full posterior saved to: {save_path}")

    def priliminary_analysis(self):

        # === MAP computation ===
        print(Fore.WHITE + "---> [INFO] Computing MAP values ...")
        # Building dictionnary
        map_values = {}
        # Gathering posteriors
        data_posteriors = self.posteriors.posterior

        # Computaing MAP
        for var in data_posteriors.data_vars:
            samples = data_posteriors[var].values.flatten()
            map_values[var] = self.get_MAP(samples)

        # === Rhat computation ===
        print(Fore.WHITE + "---> [INFO] Computing Gelman-Rubin diagnostic (R̂) ...")
        posterior = self.posteriors.posterior  
        R_hat_all = {}

        for var_name in posterior.data_vars:
            values = posterior[var_name].values  # shape: [chain, draw, *dims]

            # Dimensions : [chain, draw]
            J = values.shape[0]  # nb chains
            L = values.shape[1]  # nb draws

            # Transpose to [draw, chain] 
            X = values.T  # shape: [L, J]

            # Mean per chains
            x_j = np.mean(X, axis=0)           # [J]
            x_star = np.mean(x_j)              # scalar

            # Variance between chains
            B = (L / (J - 1)) * np.sum((x_j - x_star)**2)

            # Variance inside chains (mean of variances for each chains)
            W = np.mean(np.var(X, axis=0, ddof=1))

            # R-hat computation
            Var_hat = ((L - 1) / L) * W + (1 / L) * B
            R_hat = Var_hat / W

            R_hat_all[var_name] = R_hat
            print(Fore.MAGENTA + f"[RESULTS] R̂ for {var_name}: {R_hat:.5f}")
        
        return map_values,R_hat_all

    def get_MAP(self,samples):
        # Computing MAP
        kde = stats.gaussian_kde(samples)
        x_vals = np.linspace(min(samples), max(samples), 1000)
        densities = kde(x_vals)
        return x_vals[np.argmax(densities)]

    def plot_posteriors_custom(self, save_path,ext):
        print(Fore.BLUE + "[STEP] Plotting posterior traces & PDFs per chain (custom)")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors available. Run inference first.")
            
            if not os.path.exists(save_path):
                os.makedirs(save_path, exist_ok=True)  

            print(Fore.WHITE + "---> [INFO] Gathering posteriors data ...")

            posterior = self.posteriors.posterior

            prior_names = list(posterior.data_vars)

            n = len(prior_names)
            n_chains = posterior.sizes["chain"]
            fig, axes = plt.subplots(n, 2, figsize=(14, 2.5 * n))

            if n == 1:
                axes = [axes]

            colors = cm.tab10.colors  # palette par défaut

            print(Fore.WHITE + "---> [INFO] Plotting Posterior results ...")
            for i, name in enumerate(prior_names):
                # === Nouvelle fonctionnalité : figure individuelle ===
                fig_indiv, axs_indiv = plt.subplots(1, 2, figsize=(12, 4))
                for c in range(n_chains):
                    chain_samples = posterior[name].sel(chain=c).values.flatten()

                    axs_indiv[0].plot(chain_samples, label=f"Chain {c}", alpha=0.8, color=colors[c % len(colors)])
                    try:
                        density = stats.gaussian_kde(chain_samples)
                        x_vals = np.linspace(np.min(chain_samples), np.max(chain_samples), 300)
                        axs_indiv[1].plot(x_vals, density(x_vals), label=f"Chain {c}", color=colors[c % len(colors)])
                        axs_indiv[1].fill_between(x_vals, density(x_vals), alpha=0.2, color=colors[c % len(colors)])
                    except Exception as kde_error:
                        print(f"[WARN] KDE failed (indiv) for {name}, chain {c}: {kde_error}")

                axs_indiv[0].set_title(f"Trace - {name}")
                axs_indiv[0].set_ylabel("Value")

                axs_indiv[1].set_title(f"PDF - {name}")
                axs_indiv[1].set_ylabel("Density")

                plt.tight_layout()

                file_name = f"{name}.{ext}"
                full_path = os.path.join(save_path, file_name)
                plt.savefig(full_path, format=ext, dpi=300, bbox_inches="tight")
                print(Fore.WHITE + f"---> [INFO] Saved individual posterior plot: {file_name}")

                plt.close(fig_indiv)

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Manual plotting failed: {e}")
            sys.exit(1)

    def summary(self):
        print(Fore.BLUE + "[STEP] Generating summary of posterior distributions")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors to summarize. Run inference first.")
            summary_stats = az.summary(self.posteriors)
            print(Fore.GREEN + "---> [SUCCESS] Summary generated")
            return summary_stats
        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Summary generation failed: {e}")
            return None

    def predict_from_posterior(self, verbose=True, limit=10):
        print(Fore.BLUE + "[STEP] Predicting from posterior samples using the surrogate model")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors available. Run inference first.")
            if self.forward_model is None:
                raise RuntimeError("No surrogate model available.")

            samples = self.posteriors.posterior

            # Stacking vectors (n_features, n_samples) → (n_samples, n_features)
            stacked = []
            for name in self.model.named_vars:
                if name not in samples:
                    raise ValueError(f"Variable '{name}' not found in posterior samples.")
                stacked.append(samples[name].stack(sample=("chain", "draw")).values)

            X_samples = np.vstack(stacked).T

            # Prediction with SMT model
            y_pred = self.forward_model.predict_values(X_samples)

            print(Fore.GREEN + f"---> [SUCCESS] Surrogate model predictions completed: {len(y_pred)} samples")

            # Printing first 10 values
            if verbose:
                print(Fore.WHITE + f"--- Showing first {min(limit, len(y_pred))} predictions ---")
                for i, (x, y) in enumerate(zip(X_samples, y_pred)):
                    print(f"[#{i+1}] θ = {np.round(x, 4)} → prediction = {y[0]:.6f}")
                    if i + 1 >= limit:
                        break

            return y_pred

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Prediction from posterior failed: {e}")
            return None

# ==============================================================================================
# ------------------------
# | Bayesian model class |
# ------------------------

import pytensor.tensor as pt # type: ignore
from pytensor.graph.op import Op # type: ignore
from pytensor.graph.basic import Apply # type: ignore

# Wrapper class to use the SMT Kriging model as a custom Aesara operation
class SMStagline(Op):
    itypes = [pt.dvector]
    otypes = [pt.dscalar]

    def __init__(self, kriging_model):
        self.model = kriging_model

    def perform(self, node, inputs, outputs):
        theta = np.array(inputs[0]).reshape(1, -1)
        pred = self.model.predict_values(theta)[0][0]
        outputs[0][0] = pred                           # return scalar prediction

# Wrapper for vectorized surrogate model use (multiple predictions)
class SMStaglineVectorized(Op):
    itypes = [pt.dmatrix]  # Input: matrix of shape (n_exp, n_inputs)
    otypes = [pt.dvector]  # Output: vector of shape (n_exp,)

    def __init__(self, kriging_model):
        self.model = kriging_model

    def perform(self, node, inputs, outputs):
        X = np.array(inputs[0])  # shape = (n_exp, n_inputs)
        preds = self.model.predict_values(X).flatten()
        outputs[0][0] = preds



