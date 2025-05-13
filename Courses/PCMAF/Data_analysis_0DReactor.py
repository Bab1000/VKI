import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Loading data
target_pressure = [1500,5000,10000,20000]

for targ_p in target_pressure:

    CSV_name = f"/home/jpe/VKI/Courses/PCMAF/Results/Res_tuning_Pstat={targ_p}.csv"

    df = pd.read_csv(CSV_name)

    T = df["Tinlet"]
    dx = df["dx"]

    plt.figure()
    plt.scatter(T,np.log(dx))
    plt.xlabel("Temperature [K]")
    plt.ylabel("log(dx)")
    plt.savefig(f"/home/jpe/VKI/Courses/PCMAF/Results/Plots/Res dx_Pstat={targ_p}.jpeg",format='jpeg', dpi=300, bbox_inches='tight')





