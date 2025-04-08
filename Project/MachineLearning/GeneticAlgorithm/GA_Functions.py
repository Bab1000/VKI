import numpy as np
import sys
from colorama import Fore, Style, init
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import imageio
from tqdm import tqdm
import shutil

# Initialize colorama for colored terminal output
init(autoreset=True)



def Initialize_POP(n_p, X_Bounds, n_G=0.5, sigma_I_r=6):
    """
    Initialize Population for a genetic algorithm.

    Parameters
    ------------
    n_p : int
        Number of individuals in the population.
    X_Bounds : list
        List of tuples (min, max) defining bounds for each gene.
    n_G : float, optional (default=0.5)
        Proportion of the population initialized using a Gaussian distribution.
    sigma_I_r : float, optional (default=6)
        Interval ratio to calculate standard deviation for Gaussian initialization.
        Example: if sigma_I_r = 6, then sigma = (max - min) / 6

    Returns
    -----------
    X_V : np.ndarray (shape: number_features x n_p)
        Initial population. Each column represents an individual.
    """
    print(Fore.WHITE + "---> [INFO] Initialization of genetic algorithm")

    # Number of features/genes per individual (based on bounds)
    number_features = len(X_Bounds)

    print(Fore.WHITE + "---> [INFO] Generation of initial population")

    try:
        # ----- Gaussian-distributed portion -----
        N_Gau_pop = int(n_G * n_p)  # Number of individuals from Gaussian distribution
        X_G = np.zeros((number_features, N_Gau_pop))  # Initialize empty Gaussian population

        for j in range(number_features):
            # Compute mean and standard deviation for each gene
            mean = (X_Bounds[j][1] + X_Bounds[j][0]) / 2
            sigma = abs(X_Bounds[j][1] - X_Bounds[j][0]) / sigma_I_r
            # Generate Gaussian values for this gene across N_Gau_pop individuals
            X_G[j, :] = np.random.normal(mean, sigma, N_Gau_pop)

        # ----- Uniformly-distributed portion -----
        n_U = n_p - N_Gau_pop  # Remaining individuals from uniform distribution
        X_U = np.zeros((number_features, n_U))  # Initialize empty uniform population

        for j in range(number_features):
            # Generate uniform values between the defined min and max for each gene
            X_U[j, :] = np.random.uniform(X_Bounds[j][0], X_Bounds[j][1], n_U)

        # ----- Final assembly -----
        # Combine Gaussian and uniform subpopulations into the full population matrix
        X_V = np.concatenate([X_G, X_U], axis=1)

        print(Fore.GREEN + f"---> [SUCCESS] Initialization successfully executed")

        return X_V

    except Exception as e:
        # Handle any error that occurs during initialization
        print(Fore.RED + f"[ERROR] Error during initialization : {e}")
        sys.exit(1)  # Exit the program in case of critical error

 # 2. Evaluate Population

def Evaluate_POP(sm_q, Pdyn, Pc, Tinlet, X_V, Func, Q_target): 
    """
    Evaluate a population of candidate solutions.

    Parameters
    ------------
    sm_q : 
        Surrogate model
    Pdyn : float
        Input dynamic pressure [Pa]
    Pc : float
        Input Static pressure in the chamber [Pa]
    Tinlet : float
        Input inlet temperature [K]
    X_V : np.ndarray (shape: n_f x n_p)
        Input population. Each column is an individual (chromosome).
    Func : function
        Objective function to evaluate each individual. 
        It should take the form: Func(sm_q, Pdyn, Pc, Tinlet, individual)

    Returns
    ----------- 
    Err : np.ndarray (shape: n_p x 1)
        Cost (error) for each individual in the population.
    """
    try:
        # Extract dimensions: number of features and individuals
        n_f, n_p = X_V.shape

        # Initialize error array to store the cost for each individual
        Err = np.zeros((n_p, 1))

        # Loop over the population
        for k in range(n_p):
            # Evaluate the cost function for individual k and subtract the target
            Err[k] = abs(Func(sm_q, Pdyn, Pc, Tinlet, X_V[:, k]) - Q_target)

        return Err

    except Exception as e:
        print(Fore.RED + f"[ERROR] Population evaluation failed: {type(e).__name__} - {e}")
        sys.exit(1)

def Update_POP(X_V, Err, X_Bounds, n_I, N_ITER, mu_I=0.3, mu_F=0.05, p_M=0.5, n_E=0.05): 
    """
    Update Population based on selection, elitism, mutation, and crossover.

    Parameters
    ------------
    X_V : np.ndarray (n_f x n_p)
        Current population. Each column is an individual.
    Err : np.ndarray (n_p x 1)
        Cost (error) of each individual.
    X_Bounds : list of tuples
        Bounds (min, max) for each gene (feature).
    n_I : int
        Current generation index.
    N_ITER : int
        Total number of generations.
    mu_I : float
        Initial mutation rate.
    mu_F : float
        Final mutation rate.
    p_M : float
        Proportion of genes to mutate in each individual.
    n_E : float
        Proportion of population to keep as elites.

    Returns
    -----------
    X_V_n : np.ndarray (n_f x n_p)
        Updated population after selection and genetic operations.
    """
    try:

        # Get population shape
        n_f, n_p = X_V.shape

        # Sort individuals by fitness (ascending)
        index = Err.argsort(axis=0)

        # --- Calculate number of individuals for each genetic operation ---
        alpha = 1 / N_ITER * np.log(mu_F / mu_I)  # Exponential decay coefficient
        Mut = mu_I * np.exp(alpha * n_I)          # Current mutation rate
        N_M = int(np.round(Mut * n_p))            # Number of individuals to mutate
        N_E = int((n_p - N_M) * n_E)              # Number of elites to preserve
        N_C = int(n_p - N_M - N_E)                # Number of individuals created via crossover

        # --- 1. Elitism ---
        X_V_E = X_V[:, index[0:N_E, 0]]  # Best individuals (elites)

        # --- 2. Mutation ---
        P_M = int(p_M * n_f)             # Number of genes to mutate
        X_V_M = np.zeros((n_f, N_M))     # Container for mutated individuals

        for m in range(N_M):
            # Start with a copy of one of the best individuals
            X_V_M[:, m] = X_V[:, index[m, 0]]
            for mm in range(P_M):
                Ind_M = np.random.randint(0, n_f)  # Random gene index
                # Mutate the gene with a random value in its bounds
                X_V_M[Ind_M, m] = np.random.uniform(X_Bounds[Ind_M][0], X_Bounds[Ind_M][1])

        # --- 3. Crossover ---
        X_V_C = np.zeros((n_f, N_C))  # Container for crossover individuals

        for k in range(N_C):
            # Randomly select two parents using triangular distribution
            SEL = np.random.triangular(0, 0, N_C, 2)
            for j in range(n_f):
                a = np.random.uniform(0, 1)  # Blend factor
                # Blend two parents to create a child
                X_V_C[j, k] = a * X_V[j, index[int(SEL[0]), 0]] + (1 - a) * X_V[j, index[int(SEL[1]), 0]]

        # --- Final concatenation and boundary correction ---
        X_V_n = np.concatenate([X_V_C, X_V_E, X_V_M], axis=1)

        # Ensure all genes are within their defined bounds
        for j in range(n_f):
            mask1 = X_V_n[j, :] < X_Bounds[j][0]
            X_V_n[j, mask1] = X_Bounds[j][0]
            mask2 = X_V_n[j, :] > X_Bounds[j][1]
            X_V_n[j, mask2] = X_Bounds[j][1]

        return X_V_n

    except Exception as e:
        print(Fore.RED + f"[ERROR] Update_POP failed: {type(e).__name__} - {e}")
        return None

def Anim_COMP(sm_q, Pdyn, Pc, Tinlet, Func, X_Bounds,Q_target, n_p=100, N_ITER=100, n_G=0.5,
              sigma_I_r=6, mu_I=0.3, mu_F=0.05, p_M=0.5, n_E=0.05,
              x_1m=-2, x_1M=2, x_2m=-0.5, x_2M=3, npoints=200, Name_Video='Gif.gif'):
    """
    Animate the search process of the genetic algorithm over 2D optimization space.

    Parameters
    ----------
    - sm_q, Pdyn, Pc, Tinlet: inputs passed to the objective function `Func`
    - Func : function
        Objective function to minimize.
    - X_Bounds : list of tuples
        Bounds for each variable (min, max)
    - n_p : int
        Population size
    - N_ITER : int
        Number of generations
    - n_G, sigma_I_r, mu_I, mu_F, p_M, n_E : GA parameters
    - x_1m, x_1M, x_2m, x_2M : float
        Bounds for plotting in 2D
    - npoints : int
        Resolution of the contour plot
    - Name_Video : str
        Output filename for the resulting GIF animation

    Returns
    -------
    X_S : np.ndarray
        Best solution found
    X_U : np.ndarray
        Standard deviation of the final population
    X_V : np.ndarray
        Final population
    """
    try:
        print(Fore.BLUE + f"[STEP] Starting optimization loop")
        # Temporary folder to store frames
        FOLDER = 'Temp'
        if not os.path.exists(FOLDER):
            os.makedirs(FOLDER)

            print(Fore.WHITE + "---> [INFO] Temporary folder created")


        # Prepare the cost function surface (contour)
        x = np.linspace(x_1m, x_1M, npoints)
        y = np.linspace(x_2m, x_2M, npoints)
        X, Y = np.meshgrid(x, y)
        COST = np.zeros((npoints, npoints))

        for i in range(len(x)):
            for j in range(len(x)):
                XX = np.array([X[i, j], Y[i, j]])
                COST[i, j] = abs(Func(sm_q, Pdyn, Pc, Tinlet, XX) - Q_target)

        plt.ioff()  # Turn off interactive mode

        # Initialize population
        X_V = Initialize_POP(n_p, X_Bounds, n_G=n_G, sigma_I_r=sigma_I_r)

        # Arrays to store best and mean cost over time
        Err_Best = np.zeros((N_ITER, 1))
        Err_Mean = np.zeros((N_ITER, 1))

        for k in tqdm(range(N_ITER), desc=Fore.WHITE + f"---> [INFO] Optimization loop started"):
            # Evaluate the current population
            Err = Evaluate_POP(sm_q, Pdyn, Pc, Tinlet, X_V, Func,Q_target)

            # Update the population using genetic operators
            X_V = Update_POP(X_V, Err, X_Bounds, k, N_ITER,
                             mu_I=mu_I, mu_F=mu_F, p_M=p_M, n_E=n_E)

            # Store best and mean errors
            Err_Best[k] = np.min(np.abs(Err))
            Err_Mean[k] = np.mean(Err)

            # Plot current population on cost contour
            fig = plt.figure(figsize=(10, 4))
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)

            contour = ax1.contourf(X, Y, COST, cmap=cm.coolwarm, extend='both', alpha=0.5)
            ax1.plot(X_V[0, :], X_V[1, :], 'ko', markersize=3)
            ax1.set_xlim([x_1m, x_1M])
            ax1.set_ylim([x_2m, x_2M])
            ax1.set_xlabel(r"$\gamma_N$")
            ax1.set_ylabel(r"$\gamma_O$")                                                                                    
            # Add colorbar with label
            #cbar = plt.colorbar(contour, ax=ax1)
            #cbar.set_label("Error: prediction vs experimental", rotation=270, labelpad=15)
            


            ax2.plot(np.linspace(0, k, k), Err_Best[0:k], 'ro:', label='Best')
            ax2.plot(np.linspace(0, k, k), Err_Mean[0:k], 'bo:', label='Mean')
            ax2.legend()

            plt.title("Iteration " + str(k))
            plt.savefig(f"{FOLDER}/Step{k}.png", dpi=200)
            plt.close('all')

        # Final best result
        Index = Err.argmin()
        X_S = X_V[:, Index]
        X_U = np.std(X_V, axis=1)
        print(Fore.GREEN + f"--> [SUCCESS] Optimization successflully terminated")


        # Create the animated GIF
        print(Fore.BLUE + f"[STEP] Generation of the animation")

        images = []
        for k in tqdm(range(N_ITER), desc="Loading frames"):
            FIG_NAME = f"{FOLDER}/Step{k}.png"
            images.append(imageio.imread(FIG_NAME))

        imageio.mimsave(Name_Video, images, duration=0.5)

        # Clean up temporary folder
        shutil.rmtree(FOLDER)

        print(Fore.GREEN + f"--> [SUCCESS] Animation successfully generated")

        GammaN = X_S[0]
        GammaO = X_S[1]

        # Results
        print(Fore.MAGENTA + f"[RESULTS] GammaN = {GammaN} | GammaO = {GammaO}")

        gamma = [GammaN,GammaO]

        Q_pred = Func(sm_q, Pdyn, Pc, Tinlet, gamma)

        print(Fore.MAGENTA + f"[RESULTS] Q_pred = {Q_pred.item():.2f} | Q_exp = {Q_target}| Error = {abs(Q_pred - Q_target).item()}")

        return X_S, X_U, X_V

    except Exception as e:
        print(Fore.RED + f"[ERROR] Animation failed: {type(e).__name__} - {e}")
        sys.exit(1)

def GA_model(sm_q, Pdyn, Pc, Tinlet, Func, X_Bounds,Q_target, n_p=100, N_ITER=100, n_G=0.5,
              sigma_I_r=6, mu_I=0.3, mu_F=0.05, p_M=0.5, n_E=0.05):
    """
    Animate the search process of the genetic algorithm over 2D optimization space.

    Parameters
    ----------
    - sm_q, Pdyn, Pc, Tinlet: inputs passed to the objective function `Func`
    - Func : function
        Objective function to minimize.
    - X_Bounds : list of tuples
        Bounds for each variable (min, max)
    - n_p : int
        Population size
    - N_ITER : int
        Number of generations
    - n_G, sigma_I_r, mu_I, mu_F, p_M, n_E : GA parameters
    - x_1m, x_1M, x_2m, x_2M : float
        Bounds for plotting in 2D
    - npoints : int
        Resolution of the contour plot
    - Name_Video : str
        Output filename for the resulting GIF animation

    Returns
    -------
    X_S : np.ndarray
        Best solution found
    X_U : np.ndarray
        Standard deviation of the final population
    X_V : np.ndarray
        Final population
    """
    try:
        print(Fore.BLUE + f"[STEP] Starting optimization loop")

        # Initialize population
        X_V = Initialize_POP(n_p, X_Bounds, n_G=n_G, sigma_I_r=sigma_I_r)

        for k in tqdm(range(N_ITER), desc=Fore.WHITE + f"---> [INFO] Optimization loop started"):
            # Evaluate the current population
            Err = Evaluate_POP(sm_q, Pdyn, Pc, Tinlet, X_V, Func,Q_target)

            # Update the population using genetic operators
            X_V = Update_POP(X_V, Err, X_Bounds, k, N_ITER,
                             mu_I=mu_I, mu_F=mu_F, p_M=p_M, n_E=n_E)

        # Final best result
        Index = Err.argmin()
        X_S = X_V[:, Index]
        X_U = np.std(X_V, axis=1)
        print(Fore.GREEN + f"---> [SUCCESS] Optimization successflully terminated")

        GammaN = X_S[0]
        GammaO = X_S[1]

        # Results
        print(Fore.MAGENTA + f"[RESULTS] GammaN = {GammaN} | GammaO = {GammaO}")

        gamma = [GammaN,GammaO]

        Q_pred = Func(sm_q, Pdyn, Pc, Tinlet, gamma)

        print(Fore.MAGENTA + f"[RESULTS] Q_pred = {Q_pred.item():.2f} | Q_exp = {Q_target}| Error = {abs(Q_pred - Q_target).item()}")

        Q_pred = Func(sm_q, Pdyn, Pc, Tinlet, gamma)
        return X_S, X_U, X_V, Q_pred, abs(Q_pred - Q_target).item()

    except Exception as e:
        print(Fore.RED + f"[ERROR] Animation failed: {type(e).__name__} - {e}")
        sys.exit(1)




