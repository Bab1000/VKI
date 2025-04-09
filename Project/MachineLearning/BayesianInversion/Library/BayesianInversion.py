import pytensor  # type: ignore
pytensor.config.mode = "FAST_COMPILE"
import pymc as pm # type: ignore
import numpy as np
import arviz as az # type: ignore
import matplotlib.pyplot as plt
from colorama import Fore, Style, init
import sys
import warnings

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

            self.y = y                         # Outpu features
            self.priors = []                   # Prior distributions
            self.observed = []                 # Observed distributions
            self.model = None                  # Bayesian model
            self.forward_model = forward_model # Forward surrogate model
            self.posteriors = None             # Posterior results

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Initialization failed: {e}")
            sys.exit(1)
            raise

    def build_model(self,prior_names):
        print(Fore.WHITE + "---> [INFO] Building the Bayesian framework ...")
        try:
            with pm.Model() as model:

                # ===========================================================================================
                
                # Definition of prior distributions
                # ---------------------------------

                for dist in self.priors:
                    try:

                        # Gathering distribution infos
                        name = dist["name"]
                        dist_type = dist["type"]
                        mu = dist.get("mu", 0)
                        uncertainty = dist.get("uncertainty", 1)

                        # Construction of normal distribution
                        if dist_type == "normal":
                            sigma = max(uncertainty / 1.96, 1e-2)
                            var = pm.Normal(name, mu=mu, sigma=sigma)

                        # Construction of normal distribution
                        elif dist_type == "uniform":
                            lower = mu 
                            upper = uncertainty 
                            var = pm.Uniform(name, lower=lower, upper=upper)

                        else:
                            raise ValueError(f"Unsupported distribution: {dist_type}")

                        print(Fore.WHITE + f"---> [INFO] Added {dist_type.capitalize()} variable to model : {name}")
                    except Exception as e:
                        print(Fore.RED + f"---> [ERROR] Prior {name} distribution was not added to the model: {e}")
                        sys.exit(1)

                # ===========================================================================================

                # Forward model definition (Stagline surrogate model)
                # --------------------------------------------------

                try:
                    # Check for stagline surrogate model 
                    if self.forward_model is None:
                        raise ValueError("No surrogate model provided.")

                    # Instanciation of Stagline SM from SMStagline class to make it compatible with pyMC
                    Stagline_SM = SMStagline(self.forward_model)

                    prior_var_list = []

                    # Checking for the priors 
                    for name in prior_names:
                        # Gathering the priors just defined
                        var = model.named_vars.get(name)
                        if var is None:
                            raise ValueError(f"Prior distribution '{name}' not found in model.")
                        
                        # Gathering distributions
                        prior_var_list.append(pt.flatten(var))

                    SM_input = pt.concatenate(prior_var_list, axis=0)  # shape: (5,)
                    SM_prediction = Stagline_SM(SM_input)  # <-- dvector in, dscalar out

                    print(Fore.WHITE + "---> [INFO] Surrogate forward model integrated from 5 independent priors")

                except Exception as e:
                    print(Fore.RED + f"---> [ERROR] Failed to apply surrogate model: {e}")
                    sys.exit(1)


                # ===========================================================================================

                # Definition of observed distributions
                # ------------------------------------

                for dist in self.observed:
                    try:

                        # Gathering distribution infos
                        name = dist["name"]
                        dist_type = dist["type"]
                        mu = dist.get("mu", 0)
                        uncertainty = dist.get("uncertainty", 1)

                        # Construction of normal distribution
                        if dist_type == "normal":
                            sigma = max(uncertainty / 1.96, 1e-2)
                            # Uses Stagline SM for the observed quantities
                            var = pm.Normal(name, mu=SM_prediction, sigma=sigma, observed = self.y)

                        # Construction of normal distribution
                        elif dist_type == "uniform":
                            lower = mu 
                            upper = uncertainty 
                            var = pm.Uniform(name, lower=lower, upper=upper, observed = self.y)

                        else:
                            raise ValueError(f"Unsupported distribution: {dist_type}")

                        print(Fore.WHITE + f"---> [INFO] Added {dist_type.capitalize()} variable to model : {name}")
                    except Exception as e:
                        print(Fore.RED + f"---> [ERROR] Observed {name} distribution was not added to the model: {e}")
                        sys.exit(1)

                # Construction of the Bayesian model
                self.model = model

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Model building failed: {e}")
            sys.exit(1)
    
    def add_priors(self, name, distribution="normal", mu=None, uncertainty=None):
        try:
            distribution = distribution.lower()

            if distribution == "normal" and (mu is None or uncertainty is None):
                raise ValueError("mu and uncertainty must be provided for a Normal distribution.")
            if distribution == "uniform" and uncertainty is None:
                raise ValueError("uncertainty must be provided for a Uniform distribution.")

            # On stocke la définition de la distribution dans une liste
            self.priors.append({
                "name": name,
                "type": distribution,
                "mu": mu,
                "uncertainty": uncertainty
            })

            print(Fore.WHITE + f"---> [INFO] {distribution.capitalize()} distribution registered : {name}")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register distribution '{name}': {e}")
            sys.exit(1)

    def add_observed(self, name, distribution="normal", mu=None, uncertainty=None):
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
                "uncertainty": uncertainty
            })

            print(Fore.WHITE + f"---> [INFO] {distribution.capitalize()} distribution registered : {name}")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register distribution '{name}': {e}")
            sys.exit(1)

    def run_inference(self, prior_names ,draws = 2000, tune = 1000, chains = 4, cores = 4, restart_sim = None, save_path = None):
        print(Fore.BLUE + "[STEP] Running Bayesian inference")
        try:
            if self.model is None:
                print(Fore.WHITE + "---> [INFO] Building the model instance ...")
                self.build_model(prior_names)
                print(Fore.GREEN + "---> [SUCCESS] Model instance successfully constructed")

            # Starts the Inversion process using Markov Chain Monte Carlo process (NUTS)

            # Tune (Warm-up parameter):
            #       - MCMC has to "warm up" (adapt its internal settings like step size)
            #       - These results will not be included in the final results
            # Draws:
            #       - MCMC Sampling results

            warnings.filterwarnings("ignore", category=RuntimeWarning)
            np.seterr(over='ignore', invalid='ignore')
            
            # === Restart simulation from previous one ===
            if restart_sim:
                print(Fore.WHITE + f"---> [INFO] Resuming from saved trace: {restart_sim}")
                previous = az.from_netcdf(restart_sim)
                start = {
                    var: previous.posterior[var].mean(dim=["chain", "draw"]).values
                    for var in previous.posterior
                }

            with self.model:
                print(Fore.WHITE + "---> [INFO] Starting Bayesian inversion (MCMC sampling) ...")
                step = pm.Metropolis()
                # MCMC process
                self.posteriors = pm.sample(
                    draws=draws,                # Number of time steps per chains
                    tune=tune,                  # Number of Warm-up or burn-in steps per chains
                    chains=chains,                   # Number of chains
                    cores=cores,                    # Number of cores
                    return_inferencedata=True,
                    step=step 
                )

            print(Fore.GREEN + "---> [SUCCESS] Computation successfully completed !")

            # === Autosave ===
            if save_path:
                az.to_netcdf(self.posteriors, save_path)
                print(Fore.WHITE + f"---> [INFO] Posterior autosaved to '{save_path}'")

            return self.posteriors
        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Inference failed: {e}")

    def plot_posteriors(self):
        print(Fore.BLUE + "[STEP] Plotting the posterior samples")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors available. Run inference first.")
            az.plot_trace(self.posteriors)
            plt.savefig("Results/Posteriors.jpeg",  format='jpeg', dpi=300, bbox_inches='tight')
            print(Fore.GREEN + "---> [SUCCESS] Posteriors plot displayed !")
        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Plotting failed: {e}")
            sys.exit(1)

    def plot_posteriors_custom(self, prior_names=None, name_image = None):
        print(Fore.BLUE + "[STEP] Plotting posterior traces & PDFs per chain (custom)")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors available. Run inference first.")

            import scipy.stats as stats
            import matplotlib.cm as cm

            print(Fore.WHITE + "---> [INFO] Gathering posteriors data ...")

            posterior = self.posteriors.posterior

            if prior_names is None:
                prior_names = list(posterior.data_vars)

            n = len(prior_names)
            n_chains = posterior.sizes["chain"]
            fig, axes = plt.subplots(n, 2, figsize=(14, 2.5 * n))

            if n == 1:
                axes = [axes]

            colors = cm.tab10.colors  # palette par défaut

            print(Fore.WHITE + "---> [INFO] Plotting Posterior results ...")
            for i, name in enumerate(prior_names):
                # Accès aux échantillons par chaîne
                for c in range(n_chains):
                    chain_samples = posterior[name].sel(chain=c).values.flatten()

                    # Traceplot
                    axes[i][0].plot(chain_samples, label=f"Chain {c}", alpha=0.8, color=colors[c % len(colors)])

                    # PDF (KDE)
                    try:
                        density = stats.gaussian_kde(chain_samples)
                        x_vals = np.linspace(np.min(chain_samples), np.max(chain_samples), 300)
                        axes[i][1].plot(x_vals, density(x_vals), label=f"Chain {c}", color=colors[c % len(colors)])
                        axes[i][1].fill_between(x_vals, density(x_vals), alpha=0.2, color=colors[c % len(colors)])
                    except Exception as kde_error:
                        print(f"[WARN] KDE failed for {name}, chain {c}: {kde_error}")

                axes[i][0].set_title(f"Trace - {name}")
                axes[i][0].set_ylabel("Value")

                axes[i][1].set_title(f"PDF - {name}")
                axes[i][1].set_ylabel("Density")

            plt.tight_layout()
            plt.savefig(name_image, format="jpeg", dpi=300, bbox_inches="tight")
            print(Fore.GREEN + "---> [SUCCESS] Posterior plots Successfully saved !")

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

    def predict_from_posterior(self, prior_names, verbose=True, limit=10):
        print(Fore.BLUE + "[STEP] Predicting from posterior samples using the surrogate model")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors available. Run inference first.")
            if self.forward_model is None:
                raise RuntimeError("No surrogate model available.")

            samples = self.posteriors.posterior

            # Stacking vectors (n_features, n_samples) → (n_samples, n_features)
            stacked = []
            for name in prior_names:
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



