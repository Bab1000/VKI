import numpy as np

R = 0.03584
Pstag = 2996
rho_stag = 0.000381163

beta_f = 1/R * np.sqrt(2 * (Pstag)/rho_stag)

print(f"Beta_f = {beta_f:.2f} [1/s]")

