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
from tabulate import tabulate # type: ignore
from pytensor.printing import pydotprint # type: ignore

# Initialize colorama
init(autoreset=True)

# ==============================================================================================
# ------------------------
# | Bayesian model class |
# ------------------------

class BayesianInversion:

    def __init__(self, y, forward_model):

        """
        Initializes an instance of the Bayesian model.

        Parameters:
        - y: Observed data (target features)
        - forward_model: Surrogate model used as the forward model in Bayesian inference

        During initialization:
        - Stores the observed data
        - Prepares empty lists for prior distributions (`priors`) and observed distributions (`observed`)
        - Sets `self.model` and `self.posteriors` to None as placeholders
        - Assigns the provided surrogate model to `self.forward_model`
        """

        print(Fore.BLUE + "[STEP] Initializing the Bayesian model")
        try:
            
            print(Fore.WHITE + "---> [INFO] Building initialization ...")
        
            self.y = y                                               # Observed features
            self.priors = []                                         # Prior distributions
            self.observed = []                                       # Observed distributions
            self.model = None                                        # Bayesian model
            self.forward_model = forward_model                       # Forward surrogate model
            self.posteriors = None                                   # Posterior results
            self.SM_prediction = None                                # Link between SMT and pyMC
            self.constant = []                                       # Constant value
            self.input_names_order = ["Pdyn","Pstat","T","Gglob"]    # Input order for SM model

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Initialization failed: {e}")
            sys.exit(1)
            raise
        
    def build_model(self):

        """
        Builds the full Bayesian inference model using the provided priors, observed data, and surrogate (forward) model.

        This method performs the following steps:

        1. **Defines prior distributions**:
        Iterates over `self.priors` and adds each distribution (normal or uniform) to the PyMC model, computing standard deviations for normal priors if needed.

        2. **Integrates the surrogate model**:
        Depending on the dimensionality of `self.y`, detects whether it's a scalar (single experiment) or vector (multiple experiments).
        - For multiple experiments, groups the priors by test case and constructs a 2D input matrix to feed into the `SMStaglineVectorized` model.
        - For single experiments, concatenates the priors into a single input vector for the `SMStagline` model.
        In both cases, uses the surrogate model to predict the outcome (`SM_prediction`).

        3. **Defines observed distributions**:
        Iterates over `self.observed` to link the model prediction with the actual observed data `self.y`, using either normal or uniform distributions.

        4. **Model assignment**:
        Stores the fully constructed PyMC model in `self.model`.

        Raises:
        - ValueError for unsupported or incomplete distribution definitions
        - SystemExit if a critical step fails (e.g., prior/observed distributions not properly set or surrogate model fails)
        """


        print(Fore.WHITE + "---> [INFO] Building the Bayesian framework ...")
        try:
            with pm.Model() as model:

                # ===========================================================================================

                # Insertion of constant values
                # ----------------------------
                try:
                    if hasattr(self, "constant"):
                        for const in self.constant:
                            name = const["name"]
                            value = const["value"]
                            pm.Deterministic(name, pt.constant(value))
                except Exception as e:
                    print(Fore.RED + f"---> [ERROR] Failed to insert constants into model: {e}")
                    sys.exit(1)

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
                            sigma = max(uncertainty / 1.96)
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
                        all_row_var_names = []

                        # Gather infos for each test case
                        for tc in sorted_tc:
                            row = []
                            row_var_names = []
                            for var_name in grouped_vars[tc]:
                                # Gather each distribution of the test case
                                var = model.named_vars.get(var_name)
                                if var is None:
                                    raise ValueError(f"Prior '{var_name}' not found in model.")
                                # Appending the values in the matrix
                                vec = pt.flatten(var)[0].dimshuffle(())
                                row.append(vec)
                                row_var_names.append(var_name)  # Capture de name for further sorting



                            # Append global variables 
                            for var_name in global_vars["global"]:
                                # Gathering each global distribution 
                                var = model.named_vars.get(var_name)
                                # Appending the values in the matrix
                                vec = pt.flatten(var)[0].dimshuffle(())
                                row.append(vec)
                                row_var_names.append(var_name)  # Capture the name for further sorting

                            # Sorting the inputs for the SM model
                            index_sorted_row = []

                            for input_name in self.input_names_order:
                                found = False
                                for idx, row_var_name in enumerate(row_var_names):
                                    base_name = row_var_name.split("_TC=")[0]  # Remove TC suffix if present
                                    if input_name == base_name:
                                        index_sorted_row.append(idx)
                                        found = True
                                        break
                                if not found:
                                    raise ValueError(f"[ERROR] Input '{input_name}' not found among {row_var_names}")

                            # Now reorder row and row_var_names according to the sorted indices
                            row = [row[i] for i in index_sorted_row]
                            row_var_names = [row_var_names[i] for i in index_sorted_row]

                            # Duplicate the last element (Gglob) to satisfy inputs of the surrogate model
                            row.append(row[-1])  # duplicate the last value
                            row_var_names.append(row_var_names[-1])  # duplicate the last variable name

                            all_row_var_names.append(row_var_names)
                                                                    
                            # Creating the global input matrix for the surrogate model
                            input_matrix.append(pt.stack(row))
                        
                        # Creating a pytensor input for the SMTStaglineVectorized class
                        SM_input = pt.stack(input_matrix, axis=0)  # shape: (n_exp, n_inputs)
                        # Prediction of the output of the model
                        SM_prediction = Stagline_SM(SM_input)

                        header = [f"Input {i+1}" for i in range(len(self.input_names_order))]
                        print(Fore.WHITE + f"\n=== Input matrix for the current simulation === ")
                        print(tabulate(all_row_var_names, headers=header, tablefmt="fancy_grid"))
                        print("\n")

                    else:
                        # Non-vectorized case
                        print(Fore.WHITE + "---> [INFO] Detected scalar setup (single experiment)")
                        Stagline_SM = SMStagline(self.forward_model)

                        # Gathering the prior distributions
                        prior_var_list = []
                        prior_var_names = []  # Store names for sorting
                        for name in model.named_vars:
                            # Gathering the values for each prior distributions
                            var = model.named_vars.get(name)
                            if var is None:
                                raise ValueError(f"Prior distribution '{name}' not found in model.")
                            # Building the input vector for the SMT surrogate model
                            prior_var_list.append(pt.flatten(var))
                            prior_var_names.append(name)

                        # Sorting the inputs according to the SMT model expected order 
                        index_sorted_prior = []

                        for input_name in self.input_names_order:
                            found = False
                            for idx, prior_var_name in enumerate(prior_var_names):
                                if input_name == prior_var_name:
                                    index_sorted_prior.append(idx)
                                    found = True
                                    break
                            if not found:
                                raise ValueError(f"[ERROR] Input '{input_name}' not found among {prior_var_names}")

                        # Reorder prior_var_list
                        prior_var_list = [prior_var_list[i] for i in index_sorted_prior]
                        prior_var_names = [prior_var_names[i] for i in index_sorted_prior]

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

                    if len(np.shape(self.y)) == 1 and len(self.y) > 1: 
                        # Ensure all observed distributions are normal
                        all_normal = all(d["type"] == "normal" for d in self.observed)
                        if not all_normal:
                            raise ValueError("Vectorized implementation currently supports only normal distributions.")

                        # Compute standard deviations from 95% confidence intervals
                        sigmas = np.array([dist["uncertainty"] / 1.96 for dist in self.observed])

                        # Create a single vectorized normal distribution
                        pm.Normal("obs", mu=SM_prediction, sigma=sigmas, observed=self.y)

                        # Display mapping between symbolic prediction and observations
                        link_table = []
                        for i, dist in enumerate(self.observed):
                            name = dist["name"]
                            obs_val = self.y[i]
                            pred_expr = f"SM_prediction[{i}]"
                            link_table.append([name, str(pred_expr), obs_val / 1000])

                        print(Fore.WHITE + f"\n === Mapping between Observations and Predictions === ")
                        print(tabulate(link_table, headers=["Distribution", "SM Prediction (symbolic)", "Observed Value"], tablefmt="fancy_grid"))
                        print("\n")
                    
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
                                sigma = uncertainty / 1.96
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
                self.SM_prediction = SM_prediction
                self.model = model

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Model building failed: {e}")
            sys.exit(1)
    
    def add_priors(self, name, distribution="normal", mu=None, uncertainty=None, lower = None, upper = None):

        """
        Registers a prior distribution to be included in the Bayesian model.

        Parameters:
        - name (str): Name of the prior distribution.
        - distribution (str): Type of distribution ("normal" or "uniform"). Default is "normal".
        - mu (float, optional): Mean of the distribution (required for "normal").
        - uncertainty (float, optional): Standard deviation (required for "normal").
        - lower (float, optional): Lower bound of the distribution (required for "uniform").
        - upper (float, optional): Upper bound of the distribution (required for "uniform").

        Behavior:
        - Validates the distribution type and required parameters.
        - Appends the prior configuration to `self.priors` for later use during model construction.

        Raises:
        - ValueError: If required parameters are missing for the selected distribution type.
        """

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

            print(Fore.WHITE + f"---> [INFO] Prior distribution registered    : {name} (distribution : {distribution.capitalize()})")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register distribution '{name}': {e}")
            sys.exit(1)

    def add_observed(self, name, distribution="normal", mu=None, uncertainty=None, lower = None, upper = None):
        
        """
        Registers an observed distribution to be used in the Bayesian model.

        Parameters:
        - name (str): Name of the observed distribution.
        - distribution (str): Type of distribution ("normal" or "uniform"). Default is "normal".
        - mu (float, optional): Mean of the distribution (required for "normal").
        - uncertainty (float, optional): Standard deviation (required for "normal").
        - lower (float, optional): Lower bound (required for "uniform").
        - upper (float, optional): Upper bound (required for "uniform").

        Behavior:
        - Validates the selected distribution type and ensures required parameters are provided.
        - Appends the distribution definition to `self.observed` for later integration during model construction.

        Raises:
        - ValueError: If required parameters are missing based on the chosen distribution type.
        """

        
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

            print(Fore.WHITE + f"---> [INFO] Observed distribution registered : {name} (distribution : {distribution.capitalize()})")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register distribution '{name}': {e}")
            sys.exit(1)

    def add_constant(self, name, value):
        """
        Registers a constant value to be used in the Bayesian model.

        Parameters:
        - name (str): Name of the constant.
        - value (float): Numeric value of the constant.

        Behavior:
        - Stores the constant definition in a list for later integration.
        - The constants will be inserted into the PyMC model during `build_model()`.

        Raises:
        - ValueError: If name is missing or value is not a number.
        """

        try:
            if not isinstance(name, str):
                raise ValueError("The constant name must be a string.")
            if not isinstance(value, (int, float)):
                raise ValueError("The constant value must be a numeric type (int or float).")

            self.constant.append({
                "name": name,
                "value": float(value)
            })

            print(Fore.WHITE + f"---> [INFO] Constant registered               : {name} = {value:.3f}")

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to register constant '{name}': {e}")
            sys.exit(1)

    def run_inference(self, draws = 10000, tune = 2000, chains = 4, cores = 4, restart = False, save_path = None):

        """
        Runs the Bayesian inference process using Markov Chain Monte Carlo (MCMC) sampling.

        Parameters:
        - draws (int): Number of samples to draw from the posterior (after warm-up). Default is 10000.
        - tune (int): Number of warm-up (burn-in) steps per chain. Default is 2000.
        - chains (int): Number of independent MCMC chains to run. Default is 4.
        - cores (int): Number of CPU cores to use for parallel sampling. Default is 4.
        - restart (bool): If True, resumes sampling from a previously saved trace. Default is False.
        - save_path (str or None): Directory path to save or load the trace file (used for restarts). Default is None.

        Behavior:
        - Builds the model if it hasn't been constructed yet.
        - Optionally loads initial values from a saved trace if `restart` is enabled.
        - Calls `summary_distributions()` to display an overview of the prior and observed distributions.
        - Executes MCMC sampling using PyMC's Metropolis step method.
        - Stores the results in `self.posteriors`.
        - Automatically saves the posterior trace using `self.save_trace()` if `save_path` is provided.
        - Performs a preliminary analysis of the results with `self.priliminary_analysis()`.

        Returns:
        - posteriors (InferenceData): The full posterior trace from MCMC sampling.
        - map_values (dict): Maximum a posteriori (MAP) estimates for each parameter.
        - R_hat_all (dict): Convergence diagnostic (R-hat values) for each variable.

        Raises:
        - Exits and logs the error message if inference fails.
        """


        print(Fore.BLUE + "[STEP] Running Bayesian inference")
        try:
            if self.model is None:
                print(Fore.WHITE + "---> [INFO] Building the model instance ...")
                self.build_model()
                # Plotting symbolic tree
                #pydotprint(
                #    self.SM_prediction, 
                #    outfile="symbolic_graph.pdf", 
                #    var_with_name_simple=True, 
                #    scan_graphs=True
                #)
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
            if restart:
                
                # Check for the restart file 
                path_restart_file = os.path.join(save_path,"Restart_trace.nc")
                print(Fore.WHITE + f"---> [INFO] Resuming from saved trace: {path_restart_file}")
                previous = az.from_netcdf(path_restart_file)
                start = {
                    # Computing the mean values based on each draw of each chain
                    var: previous.posterior[var].mean(dim=["chain", "draw"]).values
                    for var in previous.posterior
                }

            # === Plotting summary of distributions ===
            self.summary_distributions()

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

    def get_constant_names(self):
        """Returns a list of all constant variable names."""
        return [const["name"] for const in self.constant]

    def priliminary_analysis(self):

        """
        Performs a preliminary analysis of the posterior results by computing:

        1. **MAP estimates**:
        - Computes the Maximum A Posteriori (MAP) value for each variable in the posterior trace.
        - Results are stored in the `map_values` dictionary.

        2. **Gelman-Rubin diagnostic (R̂)**:
        - Calculates the R-hat statistic for each variable to assess convergence across MCMC chains.
        - R̂ close to 1.0 indicates good convergence (must be < 1.2).
        - Results are stored in the `R_hat_all` dictionary and printed to the console.

        Returns:
        - map_values (dict): Dictionary containing MAP estimates for each parameter.
        - R_hat_all (dict): Dictionary of R-hat values (Gelman-Rubin diagnostic) for each parameter.

        Notes:
        - Assumes the `self.posteriors` attribute contains the posterior trace from MCMC sampling.
        - Uses custom `get_MAP()` function to extract MAP values from sample arrays.
        """


        # === MAP computation ===
        print(Fore.WHITE + "---> [INFO] Computing MAP values ...")
        # Building dictionnary
        map_values = {}
        # Gathering posteriors
        data_posteriors = self.posteriors.posterior

        # Computaing MAP
        const_names = self.get_constant_names()
        for var in data_posteriors.data_vars:
            if var in const_names:
                continue  # Skip constants
            samples = data_posteriors[var].values.flatten()
            map_values[var] = self.get_MAP(samples)

        # === Rhat computation ===
        print(Fore.WHITE + "---> [INFO] Computing Gelman-Rubin diagnostic (R̂) ...")
        posterior = self.posteriors.posterior  
        R_hat_all = {}

        for var_name in posterior.data_vars:
            if var_name in const_names:
                continue  # Skip constants

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

        """
        Estimates the Maximum A Posteriori (MAP) value from a set of posterior samples.

        Parameters:
        - samples (array-like): 1D array of posterior samples.

        Behavior:
        - Applies Gaussian Kernel Density Estimation (KDE) to approximate the probability density function of the samples.
        - Evaluates the density across a linear space of values between the min and max of the samples.
        - Returns the value corresponding to the peak (mode) of the KDE curve, which is the MAP estimate.

        Returns:
        - float: Estimated MAP value.
        """


        # Computing MAP
        kde = stats.gaussian_kde(samples)
        x_vals = np.linspace(min(samples), max(samples), 1000)
        densities = kde(x_vals)
        return x_vals[np.argmax(densities)]

    def plot_posteriors_custom(self, save_path,ext):

        """
        Generates and saves posterior plots with trace and probability density function (PDF) for each model parameter,
        including a visual overlay of the prior distribution.

        Parameters:
        - save_path (str): Directory path where the individual plots will be saved.
        - ext (str): File extension for the output figures (e.g., "png", "pdf", "svg").

        Behavior:
        - Validates that posterior samples are available; raises an error otherwise.
        - Creates the output directory if it doesn't already exist.
        - Reconstructs prior distributions from model definition (normal or uniform).
        - Iterates over each parameter in the posterior:
            - Plots trace plots of the sampled values for each chain.
            - Computes and plots the posterior density for each chain using Gaussian KDE.
            - Overlays the corresponding prior distribution as a red dashed line.
            - Automatically clips the x-axis to the support of the prior distribution:
                - μ ± 4σ for normal priors.
                - [lower, upper] bounds for uniform priors.
            - Displays a legend indicating "Posterior - Chain X" and "Prior".
        - Saves each figure to the specified path with the chosen extension.

        Raises:
        - RuntimeError: If posterior samples are not available.
        - Exits on error with detailed logging in case of plot generation or file-saving failure.
        """

        print(Fore.BLUE + "[STEP] Plotting posterior traces & PDFs per chain (custom)")
        try:
            # Ensure that posterior samples are available
            if self.posteriors is None:
                raise RuntimeError("No posteriors available. Run inference first.")

            # Create the output directory if it doesn't exist
            if not os.path.exists(save_path):
                os.makedirs(save_path, exist_ok=True)

            # Inform the user that we are loading posterior samples
            print(Fore.WHITE + "---> [INFO] Gathering posteriors data ...")

            # Load the posterior samples
            posterior = self.posteriors.posterior

            # Get the names of constant variables (excluded from plotting)
            const_names = self.get_constant_names()

            # Get all variable names except constants (i.e., parameters with priors)
            prior_names = [name for name in posterior.data_vars if name not in const_names]

            # Set color palette for different chains
            colors = cm.tab10.colors

            # Reconstruct prior distributions using scipy.stats objects
            prior_dict = {}
            for p in self.priors:
                name = p["name"]
                dist_type = p["type"]

                # Reconstruct normal prior
                if dist_type == "normal":
                    mu = p["mu"]
                    sigma = p["uncertainty"] / 1.96  # Convert 95% uncertainty to standard deviation
                    prior_dict[name] = stats.norm(loc=mu, scale=sigma)

                # Reconstruct uniform prior
                elif dist_type == "uniform":
                    lower = p["lower"]
                    upper = p["upper"]
                    prior_dict[name] = stats.uniform(loc=lower, scale=upper - lower)

            # Inform the user we are starting the plotting
            print(Fore.WHITE + "---> [INFO] Plotting Posterior results ...")

            # Get the number of MCMC chains
            n_chains = posterior.sizes["chain"]

            # Iterate through each variable
            for i, name in enumerate(prior_names):
                try:
                    # Create figure with 2 subplots (trace and PDF)
                    fig_indiv, axs_indiv = plt.subplots(1, 2, figsize=(12, 4))

                    # Store all samples from all chains for this variable
                    all_samples = []

                    # Iterate over chains to plot traces and posterior densities
                    for c in range(n_chains):
                        # Flatten the samples from chain c
                        chain_samples = posterior[name].sel(chain=c).values.flatten()
                        all_samples.append(chain_samples)

                        # Plot the trace (sample values over iterations)
                        axs_indiv[0].plot(chain_samples, label=f"Chain {c}", alpha=0.8, color=colors[c % len(colors)])

                        # Compute KDE (Kernel Density Estimation)
                        density = stats.gaussian_kde(chain_samples)

                        # Define x-axis range for KDE
                        x_vals = np.linspace(np.min(chain_samples), np.max(chain_samples), 300)

                        # Plot posterior PDF line
                        axs_indiv[1].plot(x_vals, density(x_vals), color=colors[c % len(colors)])

                        # Fill under the posterior PDF
                        axs_indiv[1].fill_between(x_vals, density(x_vals), alpha=0.2, color=colors[c % len(colors)])

                    # If a prior exists for this variable, overlay it
                    if name in prior_dict:
                        prior_dist = prior_dict[name]

                        # Identify the type of distribution
                        dist_name = getattr(prior_dist.dist, "name", "")

                        # Define support range for normal distribution
                        if dist_name == "norm":
                            x_center = prior_dist.mean()
                            x_std = prior_dist.std()
                            x_min = x_center - 4 * x_std
                            x_max = x_center + 4 * x_std

                        # Define support range for uniform distribution
                        elif dist_name == "uniform":
                            try:
                                # safer handling of loc and scale
                                loc = prior_dist.kwds.get("loc", prior_dist.args[0] if len(prior_dist.args) > 0 else -4)
                                scale = prior_dist.kwds.get("scale", prior_dist.args[1] if len(prior_dist.args) > 1 else 0)
                                x_min = loc
                                x_max = loc + scale
                            except Exception as uniform_err:
                                print(f"[WARN] Failed to extract uniform bounds for {name}: {uniform_err}")
                                # If it fails, fallback to sample-based range
                                x_min, x_max = np.min(np.concatenate(all_samples)), np.max(np.concatenate(all_samples))

                        # Fallback for any other distribution
                        else:
                            x_min, x_max = np.min(np.concatenate(all_samples)), np.max(np.concatenate(all_samples))

                        # Create x values for prior PDF
                        prior_vals = np.linspace(x_min, x_max, 500)

                        # Evaluate prior PDF
                        prior_pdf = prior_dist.pdf(prior_vals)

                        # Plot prior as a dashed red line
                        axs_indiv[1].plot(prior_vals, prior_pdf, linestyle="--", color="red", label="Prior")

                        # Add margin around the x-axis limits for clarity
                        margin = 0.02 * (x_max - x_min)
                        axs_indiv[1].set_xlim(x_min - margin, x_max + margin)

                except Exception as plot_err:
                    # Catch plotting errors for specific variables
                    print(f"[WARN] Plotting failed for {name}: {plot_err}")
                    continue

                # Set title and label for trace plot
                axs_indiv[0].set_title(f"Trace - {name}")
                axs_indiv[0].set_ylabel("Value")

                # Set title and label for PDF plot
                axs_indiv[1].set_title(f"PDF - {name}")
                axs_indiv[1].set_ylabel("Density")

                # Show legend for prior vs posterior
                axs_indiv[1].legend(loc="best")

                # Adjust layout
                plt.tight_layout()

                # Define filename and full path for saving
                file_name = f"{name}.{ext}"
                full_path = os.path.join(save_path, file_name)

                # Save figure to file
                plt.savefig(full_path, format=ext, dpi=300, bbox_inches="tight")
                print(Fore.WHITE + f"---> [INFO] Saved individual posterior plot: {file_name}")

                # Close the figure to free memory
                plt.close(fig_indiv)

        # Catch all other errors and exit safely
        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Manual plotting failed: {e}")
            sys.exit(1)

    def summary_results(self):

        """
        Generates a statistical summary of the posterior distributions using ArviZ.

        Behavior:
        - Checks if posterior data is available; raises an error if not.
        - Uses ArviZ's `summary()` function to compute key statistics (mean, std, HDI, R̂, ESS, etc.).
        - Logs success or failure with color-coded console messages.

        Returns:
        - summary_stats (DataFrame): Summary table of posterior distributions if successful.
        - None: If an error occurs during summary generation.

        Raises:
        - RuntimeError: If inference has not been run and `self.posteriors` is None.
        """


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

    def summary_distributions(self):

        """
        Displays a summary table of all prior and observed distributions defined in the model.

        Behavior:
        - Iterates through both `self.priors` and `self.observed` to extract relevant information.
        - For each distribution, displays:
        - Name
        - Role ("Prior" or "Observed")
        - Distribution type ("Normal" or "Uniform")
        - Parameters: mean, uncertainty (for normal), or lower/upper bounds (for uniform)
        - Formats and prints the summary as a well-structured table using `tabulate`.
        - Handles missing parameters by substituting placeholders ("-").

        Returns:
        - None (prints the summary directly to the console).
        """


        rows = []

        # Function to gether the format of the distribution
        def format_row(entry, role):
            # get the type of distribution
            dist_type = entry.get("type", "-").capitalize()
            # get the name
            name = entry.get("name", "-")
            
            # Normal case
            if dist_type == "Normal":
                return [
                    name,
                    role,
                    dist_type,
                    f"{entry.get('mu', '-'):.3f}" if entry.get("mu") is not None else "-",
                    f"{entry.get('uncertainty', '-'):.3f}" if entry.get("uncertainty") is not None else "-",
                    "-", "-"
                ]
            # Uniform case
            elif dist_type == "Uniform":
                return [
                    name,
                    role,
                    dist_type,
                    "-", "-",
                    f"{entry.get('lower', '-'):.3f}" if entry.get("lower") is not None else "-",
                    f"{entry.get('upper', '-'):.3f}" if entry.get("upper") is not None else "-"
                ]
            # Constant value
            elif dist_type == "Constant":
                return [
                    name,
                    role,
                    dist_type,
                    f"{entry.get('value', '-'):.3f}" if entry.get("value") is not None else "-",
                    "-", "-", "-"
                ]
            # Return something if error
            else:
                return [name, role, dist_type, "-", "-", "-", "-"]

        # Collect prior 
        for p in self.priors:
            rows.append(format_row(p, "Prior"))

        # Collect observed
        for o in self.observed:
            rows.append(format_row(o, "Observed"))

        # Constants
        for c in self.constant:
            rows.append(format_row({
                "name": c["name"],
                "type": "constant",
                "value": c["value"]
            }, "Constant"))

        # defining headers for the table
        headers = ["Name", "Role", "Distribution", "Mean", "Uncertainty", "Lower", "Upper"]

        print(Fore.WHITE + "\n=== Summary Model Distributions ===\n")
        print(Fore.WHITE + tabulate(rows, headers=headers, tablefmt="fancy_grid"))
        print("\n")

    def predict_from_posterior(self, verbose=True, limit=10):

        """
        Generates predictions using the surrogate model based on samples from the posterior distribution.

        Parameters:
        - verbose (bool): If True, prints the first few predictions to the console. Default is True.
        - limit (int): Maximum number of predictions to print when verbose is enabled. Default is 10.

        Behavior:
        - Ensures posterior samples and the surrogate model are available.
        - Stacks posterior samples into a 2D array with shape (n_samples, n_features).
        - Uses the surrogate model (`self.forward_model`) to predict output values from the sampled parameters.
        - Optionally prints a limited number of predictions for inspection.

        Returns:
        - y_pred (ndarray): Predicted outputs for each posterior sample (shape: [n_samples, 1]).

        Raises:
        - RuntimeError: If posterior samples or the surrogate model are not available.
        - ValueError: If a required variable is missing in the posterior sample set.
        - Returns None in case of failure with an error message.
        """

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

    def plot_predictions_vs_observed(self, n_samples=100, save_path=None):
        """
        Trace les prédictions du modèle SMT à partir des échantillons postérieurs
        et les compare aux observations réelles.
        
        Parameters:
        - n_samples (int): nombre d'échantillons postérieurs à utiliser pour la prédiction
        - save_path (str or None): si fourni, sauvegarde la figure dans ce chemin
        """
        print(Fore.BLUE + "[STEP] Plotting SMT predictions vs observed data")
        try:
            if self.posteriors is None:
                raise RuntimeError("No posteriors available.")
            if self.forward_model is None:
                raise RuntimeError("No surrogate model available.")

            samples = self.posteriors.posterior
            prior_names = list(self.model.named_vars)

            # === Extraction d'échantillons postérieurs
            stacked = []
            for name in prior_names:
                val = samples[name].stack(sample=("chain", "draw")).values
                stacked.append(val)
            X_all = np.vstack(stacked).T  # shape (n_total_samples, n_inputs)

            # === Réduction à n_samples max
            idx = np.random.choice(X_all.shape[0], size=min(n_samples, X_all.shape[0]), replace=False)
            X_samples = X_all[idx, :]

            # === Prédiction avec le modèle SMT
            y_preds = self.forward_model.predict_values(X_samples)  # shape: (n_samples, n_outputs)

            # === Tracé pour chaque sortie
            n_outputs = y_preds.shape[1] if y_preds.ndim == 2 else 1
            y_obs = np.atleast_2d(self.y)  # Ensure shape (n_exp, n_outputs)

            if n_outputs == 1:
                y_preds = y_preds.reshape(-1, 1)

            fig, axes = plt.subplots(n_outputs, 1, figsize=(8, 4 * n_outputs))
            if n_outputs == 1:
                axes = [axes]

            for j in range(n_outputs):
                ax = axes[j]
                preds_j = y_preds[:, j]  # shape: (n_samples,)
                obs_j = y_obs[:, j]

                pred_mean = preds_j.mean()
                pred_hdi = np.percentile(preds_j, [2.5, 97.5])

                ax.axhline(y=obs_j[0], color="black", linestyle="--", label="Observed")
                ax.hist(preds_j, bins=30, alpha=0.6, color="skyblue", label="Posterior preds")
                ax.axvline(pred_mean, color="blue", linestyle="-", label="Mean prediction")
                ax.axvspan(pred_hdi[0], pred_hdi[1], alpha=0.2, color="blue", label="95% HDI")

                ax.set_title(f"Prediction vs Observed (Output #{j+1})")
                ax.set_xlabel("Predicted value")
                ax.set_ylabel("Frequency")
                ax.legend()

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                print(Fore.WHITE + f"---> [INFO] Plot saved to: {save_path}")
            else:
                plt.show()

        except Exception as e:
            print(Fore.RED + f"---> [ERROR] Failed to plot prediction vs observed: {e}")


# ==============================================================================================
# ------------------------
# | Bayesian model class |
# ------------------------

import pytensor.tensor as pt # type: ignore
from pytensor.graph.op import Op # type: ignore
from pytensor.graph.basic import Apply # type: ignore

# Wrapper class to use the SMT Kriging model as a custom Aesara operation
class SMStagline(Op):

    """
    PyTensor-compatible wrapper for a scalar-output surrogate (Kriging) model.

    This custom `Op` allows a surrogate model to be integrated into a PyMC probabilistic graph.

    Attributes:
    - itypes: Input types for the Op (`pt.dvector`) : expects a 1D input vector of parameters.
    - otypes: Output type of the Op (`pt.dscalar`) : returns a single scalar prediction.

    Parameters:
    - kriging_model: A trained surrogate model (e.g., Kriging) with a `predict_values()` method.

    Methods:
    - perform(): Reshapes the input vector into a 2D array, runs the surrogate model prediction,
    and returns the scalar output to be used inside the PyMC model.
    """


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

    """
    PyTensor-compatible wrapper for a vector-output surrogate (Kriging) model handling multiple experiments.

    This custom `Op` allows batch predictions from a surrogate model to be integrated into a PyMC model for
    vectorized inference.

    Attributes:
    - itypes: Input types for the Op (`pt.dmatrix`) : expects a 2D input matrix of shape (n_exp, n_inputs),
    where each row corresponds to a different experiment.
    - otypes: Output type of the Op (`pt.dvector`) : returns a 1D array of scalar predictions (one per experiment).

    Parameters:
    - kriging_model: A trained surrogate model (e.g., Kriging) with a `predict_values()` method that supports batch input.

    Methods:
    - perform(): Receives a matrix of input samples, performs batch predictions using the surrogate model,
    flattens the output, and returns a vector of predicted values.
    """


    itypes = [pt.dmatrix]  # Input: matrix of shape (n_exp, n_inputs)
    otypes = [pt.dvector]  # Output: vector of shape (n_exp,)

    def __init__(self, kriging_model):
        self.model = kriging_model

    def perform(self, node, inputs, outputs):
        X = np.array(inputs[0])  # shape = (n_exp, n_inputs)
        preds = self.model.predict_values(X).flatten()
        outputs[0][0] = preds



