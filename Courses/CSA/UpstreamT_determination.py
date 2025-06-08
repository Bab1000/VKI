import numpy as np
from scipy.optimize import root_scalar
from Utils_UpstreamT_determination import UpstreamT
 
# Constants
R_univ = 8.314462618 # J/mol/K
H0_target = 18.59e6  # J/kg
M1 = 8.231
gamma = 1.314
Pstat = 74.51

# List of species for air_5
list_species = ["O", "O2", "N", "NO", "N2"]

T_sol,R = UpstreamT(M1,H0_target,list_species,gamma,Pstat)

V1 = M1 * np.sqrt(gamma*R*T_sol)

print(f"V1 = {V1:.2f}")