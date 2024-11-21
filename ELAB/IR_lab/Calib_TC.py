import numpy as np
#import fig_management
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.style.use("grayscale")

def func(x,a,b):
    return a*x + b
def func2(x,a,b):
     return (x-b)/a

# Data calib for report (if these data are commented dont forget to comment the figure saving in the plots)
T = np.array([23.4, 32.6, 37.4, 41.6, 46.1, 50.4, 55.2, 59.9, 64.4, 68.8, 59, 49.9, 40.6])
TC1 = np.array([26.17, 35.13, 40.01, 44.43, 48.78, 53.02, 57.87, 62.45, 67.04, 71.39, 62.04, 52.38, 43.19])
TC2 = np.array([26.115, 35.17, 40.03, 44.45, 48.8, 53.07, 57.91, 62.5, 67.1, 71.46, 62.07, 52.39, 43.21])

# Data calib for computation
# T = np.array([31.8, 42.2, 50.3, 59.0, 68.6])
# TC1 = np.array([36.32, 46.16, 53.26, 63.38, 72.43])
# TC2 = np.array([36.37, 46.23, 53.33, 63.45, 72.5])

T_err = np.ones(np.size(T))*0.1

coef1, err1 = curve_fit(func, TC1, T, p0=[1,0],sigma=T_err)
coef2, err2 = curve_fit(func, TC2, T, p0=[1,0],sigma=T_err)
err1 = np.sqrt(np.diag(err1))
err2 = np.sqrt(np.diag(err2))

print(f"TC1 coefficients: a = {coef1[0]:.4f} ± {err1[0]:.4f}, b = {coef1[1]:.4f} ± {err1[1]:.4f}")
print(f"TC2 coefficients: a = {coef2[0]:.4f} ± {err2[0]:.4f}, b = {coef2[1]:.4f} ± {err2[1]:.4f}")


T_dummy = np.linspace(20, 75, 1000)

plt.figure(figsize=(8, 6), dpi=100)
plt.plot(TC1[:10], T[:10], "o", fillstyle="none", ms=8, alpha=0.8, label="Increasing $T$")
plt.plot(TC1[10:], T[10:], "v", fillstyle="none", ms=8, alpha=0.8, label="Decreasing $T$")
plt.plot(T_dummy, func(T_dummy, *coef1), color="black", linewidth=1.5, linestyle="--", label="Linear Regression")
plt.legend(fontsize=10, loc="upper left", frameon=False)
plt.xlabel("$T_{\mathrm{C}_1}$ [$^{\\circ}$C]", fontsize=12)
plt.ylabel("$T$ [$^{\\circ}$C]", fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(color="gray", linestyle=":", linewidth=0.5, alpha=0.7)
plt.savefig("/home/jpe/VKI/ELAB/IR_lab/Plots/TC1_calibration_curve.pdf", bbox_inches="tight", dpi=300)
#plt.show()

plt.figure(figsize=(8, 6), dpi=100)
plt.plot(TC2[:10], T[:10], "o", fillstyle="none", ms=8, alpha=0.8, label="Increasing $T$")
plt.plot(TC2[10:], T[10:], "v", fillstyle="none", ms=8, alpha=0.8, label="Decreasing $T$")
plt.plot(T_dummy, func(T_dummy, *coef1), color="black", linewidth=1.5, linestyle="--", label="Linear Regression")
plt.legend(fontsize=10, loc="upper left", frameon=False)
plt.xlabel("$T_{\mathrm{C}_2}$ [$^{\\circ}$C]", fontsize=12)
plt.ylabel("$T$ [$^{\\circ}$C]", fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(color="gray", linestyle=":", linewidth=0.5, alpha=0.7)
plt.savefig("/home/jpe/VKI/ELAB/IR_lab/Plots/TC2_calibration_curve.pdf", bbox_inches="tight", dpi=300)