import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data_path = "Air_5_SEB_BIG_CSA_flowfield.dat"
data = np.loadtxt(data_path, skiprows=2)  # Skip the first two rows as they are headers

# Extract the columns
xc = data[:, 0]
T = data[:, 3]
p = data[:, 4]
rho = data[:, 5]



# Plot xc vs T
plt.figure(figsize=(10, 6))
plt.plot(xc, T, label='Temperature (T)', color='r')
plt.xlabel('xc')
plt.ylabel('T')
plt.title('Plot of xc vs. T')
plt.grid(True)
plt.legend()
plt.show()

# Plot xc vs p
plt.figure(figsize=(10, 6))
plt.plot(xc, p, label='Pressure (p)', color='b')
plt.xlabel('xc')
plt.ylabel('p')
plt.title('Plot of xc vs. p')
plt.grid(True)
plt.legend()
plt.show()

# Plot xc vs rho
plt.figure(figsize=(10, 6))
plt.plot(xc, rho, label='Density (rho)', color='g')
plt.xlabel('xc')
plt.ylabel('rho')
plt.title('Plot of xc vs. rho')
plt.grid(True)
plt.legend()
plt.show()


