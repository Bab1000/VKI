# %%
import numpy as np
#import fig_management
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def func(x,a,b):
    return a*x + b
def func2(x,a,b):
     return (x-b)/a

# %%
T = np.array([65.1, 60.9, 55.9, 50.3, 44.5, 40, 21.8, 26.2, 30.3, 35.5, 40.2, 44.4, 51, 54.8, 60.1, 63.8, 68.7])
TC1 = np.array([65, 60.4, 56.2, 50.3, 44.5, 40.1, 21.8, 26.4, 30.6, 35.1, 40.1, 44.4, 50.7, 54.7, 59.8, 63.5, 68.4])
TC2 = np.array([65, 60.6, 55.7, 50.6, 44.8, 40.3, 21.7, 26.3, 30.4, 35, 40.5, 44.6, 50.8, 54.7, 59.6, 63.4, 68.6])
T_err = np.ones(np.size(T))*0.1

coef1, err1 = curve_fit(func, TC1, T, p0=[1,0],sigma=T_err)
coef2, err2 = curve_fit(func, TC2, T, p0=[1,0],sigma=T_err)
err1 = np.sqrt(np.diag(err1))
err2 = np.sqrt(np.diag(err2))

print(coef1, err1)
print(coef2, err2)

#%%
plt.figure()

T_dummy = np.linspace(20,70,1000)
plt.plot(TC1[6:],T[6:],"k^",fillstyle="none",ms=5,label="Increasing $T$")
plt.plot(TC1[:6],T[:6],"kv",fillstyle="none",ms=5,label="Decreasing $T$")
plt.plot(T_dummy,func(T_dummy,*coef1),label="Linear Regression")
plt.legend()
plt.xlabel("$T_{\mathrm{C}_1}$ [$^{\\circ}$C]")
plt.ylabel("$T$ [$^{\\circ}$C]")
plt.savefig("TC1_calibration_curve.pdf")

plt.figure()
plt.plot(TC2[6:],T[6:],"k^",fillstyle="none",ms=5, label="Increasing $T$")
plt.plot(TC2[:6],T[:6],"kv",fillstyle="none",ms=5, label="Decreasing $T$")
plt.plot(T_dummy,func(T_dummy,*coef2), label="Linear Regression")
plt.xlabel("$T_{\mathrm{C}_2}$ [$^{\\circ}$C]")
plt.ylabel("$T$ [$^{\\circ}$C]")
plt.legend()
plt.savefig("TC2_calibration_curve.pdf")
# %%
