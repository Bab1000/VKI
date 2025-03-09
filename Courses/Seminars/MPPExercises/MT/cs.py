import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({'font.size': 14})

# Opening file
#US3D 


dat_i67 = np.loadtxt('./MPP.dat' , skiprows = 1)




#MPP
points = len(dat_i67)
#points = 800
time_mpp = dat_i67[0:points, 0]
csO2_mpp = dat_i67[0:points, 1]
csO_mpp  = dat_i67[0:points, 2]

# Y
plt.figure(figsize=(6,5))
plt.plot(time_mpp, csO2_mpp, linestyle='-', color='orange', linewidth=1, label=r'$Y(O2)$ MPP ODE non pref')#, label=r'$Y_{O2}$ mpp')
plt.plot(time_mpp, csO_mpp, linestyle='-', color='green' , linewidth=1, label=r'$Y(O) $ MPP ODE non pref')#, label=r'$Y_{O}$ mpp')
plt.xscale("log")
#plt.xticks([0, 4.E-6, 8.E-6])
#plt.yticks([0, 0.4, 0.8])
plt.xlabel('time [s]')
plt.ylabel('Y')
plt.legend(loc=0)
plt.savefig("Y.png", format="png", bbox_inches="tight")
plt.show()

