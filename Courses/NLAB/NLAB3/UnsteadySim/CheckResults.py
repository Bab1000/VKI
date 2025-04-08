import pandas as pd
import matplotlib.pyplot as plt
 
# Read the data
file_path = "history.csv"  # Change this to your file path
df = pd.read_csv(file_path, skipinitialspace=True)
df.columns = df.columns.str.strip()
 
# Plotting
plt.figure(figsize=(10, 6))
plt.plot(df['rms[P]'])

plt.xlabel('Inner_Iter')
plt.ylabel('Values')
plt.title('Data Visualization')
plt.legend()
plt.grid(True)
#plt.ylim(0,4)
plt.show()
 