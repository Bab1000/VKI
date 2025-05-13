import numpy as np
from scipy.optimize import root_scalar
from Utils_UpstreamT_determination import UpstreamT
 
# Constants
R_univ = 8.314462618 # J/mol/K
H0_target = 18.59e6  # J/kg
M1 = 6
gamma = 1.3237

# List of species for air_5
list_species = ["O", "O2", "N", "NO", "N2"]

T_sol = UpstreamT(M1,H0_target,list_species,gamma)