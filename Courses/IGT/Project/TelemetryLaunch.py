import numpy as np
import pandas as pd
from nrlmsise00 import msise_model
from datetime import datetime, timedelta
import matplotlib.pyplot as plt 

# ===========================================================================================

# ----------------
# | MISSION NAME |
# ----------------

mission_name = "SES-11"
launch_date = datetime(2017, 5, 1)

# ===========================================================================================

# ----------------
# | LOADING DATA |
# ----------------

csv_file_name = mission_name + "Stage1raw.xlsx"

df = pd.read_excel(csv_file_name)

time = np.array(df["time"]) # [s]
velocity = np.array(df["velocity"]) # [m/s]
altitude = np.array(df["altitude"]) # [km]

# ===========================================================================================

# --------------------------
# | ATMOSPHERIC PROPERTIES |
# --------------------------
density_list = []
a_list = []

f107 = 150 # valeur typique moyenne(intensité solaire)
f107A = 150    # moyenne 81 jours
ap = 4 # calme géomagnétique

for alt, t in zip(altitude, time):
    current_time = launch_date + timedelta(seconds=float(t))
    density,temperature = msise_model(current_time, alt, 0, 0, f107, f107A, ap)
    temp = temperature[1]  # Kelvin
    density_list.append(density[5])  # kg/m³
    a_list.append(np.sqrt(1.4 * 287 * temp))  # m/s

Mach = velocity/a_list
#Re = density_list*velocity*3.7/mu

plt.figure()
plt.plot(time,velocity)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
plt.savefig(f"/home/jpe/VKI/Courses/IGT/Project/Plots/VelvsT{mission_name}.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)

#plt.figure()
#plt.plot(time,density_list)
#plt.show()
#
plt.figure()
plt.plot(time,Mach)
plt.xlabel("Time [s]")
plt.ylabel("Mach [-]")
plt.savefig(f"/home/jpe/VKI/Courses/IGT/Project/Plots/MachvsT{mission_name}.jpeg", format='jpeg', dpi=600, bbox_inches='tight', transparent=True)













