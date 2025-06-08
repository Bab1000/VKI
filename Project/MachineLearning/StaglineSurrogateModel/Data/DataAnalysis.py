import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Path of Training data
XTrainingData_path = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Data/Training_Xdata.dat" 
YTrainingData_path = "/home/jpe/VKI/Project/MachineLearning/StaglineSurrogateModel/Data/YTraining.dat" 

# Reading inputs training data 
df = pd.read_csv(XTrainingData_path, sep=r'\s+')
Temperatures = df["Tedg[K]"].values
gammaN = df["GN"]
gammaO = df["GO"]
Pstat= df["psta[Pa]"]
Pdyn= df["pdyn[Pa]"]

print(f"Max GN training = {max(gammaN)}")
print(f"Max GO training = {max(gammaO)}")
print(f"Max Pstat training = {max(Pstat)}")
print(f"Min Pstat training = {min(Pstat)}")
print(f"Max T training = {max(Temperatures)}")
print(f"Min T training = {min(Temperatures)}")
print(f"Max Pdyn training = {max(Pdyn)}")
print(f"Min Pdyn training = {min(Pdyn)}")


# Reading outputs training data 
df = pd.read_csv(YTrainingData_path, sep=r'\s+')
qwall = df["qwall"].values/1000
residuals = df["res"].values

# Concatenating both vectors
Res_vec = np.column_stack((Temperatures, qwall,residuals,Pstat))

# Mask for NotConv values
mask = Res_vec[:,2] != "NotConv"
Res_vec = Res_vec[mask]

# Sort values of temperature
Res_vec_sorted = Res_vec[Res_vec[:, 0].argsort()]

# Dropping non-converged values
mask = Res_vec_sorted[:,2].astype(float) <= -2
Res_vec_sorted = Res_vec_sorted[mask]

print(f"Total training points:{len(Res_vec_sorted)}")

# ===========================================

plt.figure()

mask_Pstat = (Res_vec_sorted[:,3].astype(float) >= 500) & (Res_vec_sorted[:,3].astype(float) <= 2000)

# Affichage des points filtrés
plt.scatter(Res_vec_sorted[mask_Pstat, 3], Res_vec_sorted[mask_Pstat, 0], s=10, c='blue', alpha=0.7)

plt.xlabel("Static pressure [Pa]")
plt.ylabel("Temperature [K]")
plt.savefig("PstatvsT.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
plt.close()



# ===========================================


# Plot the results
plt.figure()
plt.scatter(Res_vec_sorted[:,0],Res_vec_sorted[:,1],s=10, c='blue', alpha=0.7)
plt.xlabel("Temperature [K]")
plt.ylabel("Heat flux [kW]")
plt.savefig("TvsQ.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
plt.close()

# Plot the results
plt.figure()
plt.scatter(np.log10(gammaN),np.log10(gammaO),s=10, c='blue', alpha=0.7)
plt.xlabel("log10(Gamma N)")
plt.ylabel("log10(Gamma O)")
plt.savefig("GNvsGO.jpeg", format='jpeg', dpi=300, bbox_inches='tight')
plt.close()



