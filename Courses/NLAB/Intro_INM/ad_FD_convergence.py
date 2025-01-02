import numpy as np
import matplotlib.pyplot as plt
import time



from ad_FD_1 import advect_gaussian_vectorized
from ad_FD_1 import gaussian
from ad_FD_1 import advect_gaussian_vectorized_second_order
from ad_FD_1 import plot_evolution



n_elem=[32,64,128,256,512]
e_1=np.zeros((len(n_elem),1))
e_2=np.zeros((len(n_elem),1))
errors_1st=[]
errors_2nd=[]
mesh_sizes=[]
for i in range(len(n_elem)):

    nx=n_elem[i]+1
    # Simulation parameters
    L = 1.0        # Length of the domain
    # Definition of the time step based on CFL
    CFL=2.0
    c = 1.0        # Advection speed
    tfinal=1e-1

    # Spatial grid
    dx = L / (nx - 1)
    print("--- Compute advection for ",n_elem[i], " with dx = ",dx," the dt is using CFL= ",CFL)
    x = np.linspace(0, L, nx)

    # Initial condition: Gaussian centered at x0
    sigma = 0.1 * L
    x0 = L*0.5  # Initial position of the Gaussian peak
    u = gaussian(x, x0, sigma)


    dt=CFL*dx/(c)
    nt=int(tfinal/dt)

    u_vect=advect_gaussian_vectorized(u,x,dt,nt,c)

    u_ref = gaussian(x, x0+c*nt*dt, sigma)

    u_2nd=advect_gaussian_vectorized_second_order(u,x,dt,nt,c)

    #print('Stable time step for CFL=',CFL,' is dt=',dt, ' number of time steps nt=',nt)
    #plot_evolution(x, u_2nd, nt, dt)
    #plt.figure(figsize=(10, 6))
    #plt.plot(x, np.array(u_vect)[-1,:],label="solution")
    #plt.plot(x, u,label="initail")
    #plt.plot(x, u_ref,label="reference")
    #plt.legend()
    #plt.grid(True)

    #Computing the integral of the error
    error_1st = np.linalg.norm(np.array(u_vect)[-1,:] - u_ref, ord=2) * dx
    error_2nd = np.linalg.norm(np.array(u_2nd)[-1,:] - u_ref, ord=2) * dx

    mesh_sizes.append(dx)
    errors_1st.append(error_1st)
    errors_2nd.append(error_2nd)


plt.figure(figsize=(10, 6))
plt.loglog(mesh_sizes, errors_1st, 'o-', label='Error upwind 1st')
plt.loglog(mesh_sizes, errors_2nd, 'o-', label='Error upwind 2nd')
plt.xlabel('Mesh size (h)')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
convergence_rate = np.polyfit(np.log(mesh_sizes), np.log(errors_1st), 1)[0]
print(f"Convergence rate for upwind 1st scheme: {convergence_rate}")
convergence_rate = np.polyfit(np.log(mesh_sizes), np.log(errors_2nd), 1)[0]
print(f"Convergence rate for upwind 2nd scheme: {convergence_rate}")

plt.show()

