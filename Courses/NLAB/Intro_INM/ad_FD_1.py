import numpy as np
import matplotlib.pyplot as plt
import time
import argparse

# Function to initialize a Gaussian profile
def gaussian(x, x0, sigma):
    return np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))

# Function to apply finite difference for the advection equation
def advect_gaussian(u0,mesh,dt,nt,c):
    dx=mesh[1]-mesh[0]
    # Prepare to store solutions at each time step
    u_all = [u0.copy()]

    # Time stepping loop (Euler explicit)
    for n in range(nt):
        u_new = u0.copy()

        # Finite difference scheme (upwind)
        for i in range(1, len(mesh)):
            u_new[i] = u0[i] - c * dt / dx * (u0[i] - u0[i - 1])

        # Periodic boundary condition
        u_new[0] = u0[0] - c * dt / dx * (u0[0] - u0[-1])

        u0 = u_new.copy()
        u_all.append(u0.copy())

    return u_all

def advect_gaussian_vectorized(u0,mesh,dt,nt,c):

    # Prepare to store solutions at each time step
    u_all = [u0.copy()]
    dx=mesh[1]-mesh[0]

    # Precompute constants to avoid repeating calculations
    coeff = c * dt / dx

    # Time stepping loop (Euler explicit), avoiding internal loops
    for n in range(nt):
        # Apply periodic boundary condition using NumPy roll
        u_new = u0 - coeff * (u0 - np.roll(u0, 1))

        u0 = u_new.copy()
        u_all.append(u0)

    return u_all


def advect_gaussian_vectorized_second_order(u0, mesh, dt, nt, c):
    dx = mesh[1] - mesh[0]  # Assume uniform spacing
    u_all = [u0.copy()]  # To store the solution at each time step

    # Compute the Courant number
    coeff = c * dt / dx
    t=0
    # Time-stepping loop
    for n in range(nt):
        # Using the second-order upwind scheme
        u_new = u0 - coeff * (1.5 * u0 - 2.0 * np.roll(u0, 1) + 0.5 * np.roll(u0, 2))

        # Update for the next iteration
        u0 = u_new.copy()
        u_all.append(u0)
        t+=dt
    return u_all


def advect_gaussian_vectorized_initialized(u0,mesh,dt,nt,c):
    dx = mesh[1] - mesh[0]  # Assume uniform spacing
    # Prepare to store solutions at each time step
    u_all = np.zeros((len(u0),nt+1))
    u_all[:,0]=u0

    # Precompute constants to avoid repeating calculations
    coeff = c * dt / dx

    # Time stepping loop (Euler explicit), avoiding internal loops
    for n in range(1,nt+1):
        # Apply periodic boundary condition using NumPy roll
        u_new = u0 - coeff * (u0 - np.roll(u0, 1))

        u0 = u_new.copy()
        u_all[:,n]=u0

    return u_all



def computeL2Norm(u,u_ref):
    error=np.sqrt(np.sum((u - u_ref) ** 2))
    return error

# Function to plot the evolution of the Gaussian profile
def plot_evolution(x, u_all, nt, dt):
    plt.figure(figsize=(10, 6))
    for i, u in enumerate(u_all[::int(nt / 10)]):  # plot every 10% of the time steps
        plt.plot(x, u, label=f't={i*dt*nt/10:.2f}s')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.title('Evolution of the Gaussian profile')
    plt.legend()
    plt.grid(True)

# Define the source term using analytical solution
def evaluate_MMS_source():
    source_term=0.0
    return source_term

def main(nx,tfinal,CFL,case):
    # Start measuring CPU time

    # Simulation parameters
    L = 1.0        # Length of the domain
    

    # Definition of the time step based on CFL
    c = 1.0        # Advection speed

    # Spatial grid
    dx = L / (nx - 1)
    x = np.linspace(0, L, nx)

    # Initial condition: Gaussian centered at x0
    sigma = 0.1 * L
    x0 = L*0.5  # Initial position of the Gaussian peak
    u = gaussian(x, x0, sigma)


    dt=CFL*dx/(c)
    nt=int(tfinal/dt)

    print('Stable time step for CFL=',CFL,' is dt=',dt, ' number of time steps nt=',nt)

    start_cpu_time = time.process_time()
    # Run the simulation
    if(case==1):
        u_loop = advect_gaussian(u,x,dt,nt,c)

    # Stop measuring CPU time
    end_cpu_time_loop = time.process_time()
    if(case==1):
        cpu_time_spent_loop = end_cpu_time_loop - start_cpu_time
        print(f"CPU time spent using loops: {cpu_time_spent_loop} seconds")
    if(case==2):

        u_vect=advect_gaussian_vectorized(u,x,dt,nt,c)
        end_cpu_time_vector = time.process_time()
        cpu_time_spent_vector = end_cpu_time_vector - end_cpu_time_loop
        print(f"CPU time spent using vectorized: {cpu_time_spent_vector} seconds")

        # u_vect_mem=advect_gaussian_vectorized_initialized(u,x,dt,nt,c)
        # end_cpu_time_vector_mem = time.process_time()
        # cpu_time_spent_vector_mem = end_cpu_time_vector_mem - end_cpu_time_vector
        # print(f"CPU time spent using vectorized and initialization of the memory: {cpu_time_spent_vector_mem} seconds")


    # Plot the results
    if(case==1):
        plot_evolution(x, u_loop, nt, dt)
    if(case==2):
        plot_evolution(x, u_vect, nt, dt)

    if(case==3):
        u_2nd=advect_gaussian_vectorized_second_order(u,x,dt,nt,c)
        plot_evolution(x, u_2nd, nt, dt)
    
    plt.show()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=float,help='final time', required=False, default=1.0)
    parser.add_argument('-CFL', type=float, help='CFL', required=False, default=0.1)
    parser.add_argument('-nx', type=int, help='Number of nodes (elem+1)', required=False, default=1000)
    parser.add_argument('-case', type=int, help='1 for loop, 2 for comparison CPU, 3 for second order upwind', required=False, default=3)
    args = parser.parse_args()
    nx=args.nx
    CFL=args.CFL
    tfinal=args.t
    case = args.case
    main(nx,tfinal,CFL,case)


