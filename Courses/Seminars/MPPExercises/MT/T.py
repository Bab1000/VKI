import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

#MPP ODE
dat_ODE = np.loadtxt('./MPP.dat' , skiprows = 1)
points   = len(dat_ODE)
time_mpp_ODE = dat_ODE[0:points, 0]
T_mpp_ODE    = dat_ODE[0:points, 3]
Tv_mpp_ODE   = dat_ODE[0:points, 4]


# Y
plt.figure(figsize=(6,5))
plt.plot(time_mpp_ODE, T_mpp_ODE , linestyle='-' , color='black' , linewidth=1, label='TT, MPP ODE')
plt.plot(time_mpp_ODE, Tv_mpp_ODE, linestyle='-' , color='orange', linewidth=1, label='TV, MPP ODE')
plt.xscale("log")
#plt.xticks([0, 4.E-6, 8.E-6])
#plt.yticks([0, 0.4, 0.8])
plt.xlabel('time [s]')
plt.ylabel('T [K]')
plt.legend(loc=0)
plt.savefig("T.png", format="png", bbox_inches="tight")
plt.show()

