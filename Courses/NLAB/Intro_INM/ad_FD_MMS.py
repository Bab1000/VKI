# Method of manufacturing solution
from sympy import symbols, diff, sin, cos
import numpy as np
import matplotlib.pyplot as plt

def compute_manufactured_solution(u,c, x_var, t_var):
    """
    Compute symbolic derivatives of an expression with respect to x and t.

    Parameters:
    - expression: sympy expression to differentiate.
    - x_var: sympy symbol for the variable x.
    - t_var: sympy symbol for the variable t.
    - order: integer specifying the order of derivative (default is 1 for first derivative).

    Returns:
    - A dictionary with derivatives with respect to x and t up to the specified order.
    """
    source_term=diff(u, t_var, 1)+c*diff(u, x_var, 1)
    print("Source term :",source_term)
    return source_term


def advect_vectorized(u0,mesh,dt,nt,c,source_term):

    # Prepare to store solutions at each time step
    u_all = [u0.copy()]
    dx=mesh[1]-mesh[0]

    # Precompute constants to avoid repeating calculations
    coeff = c * dt / dx

    # Time stepping loop (Euler explicit), avoiding internal loops
    ti=0.0;
    for n in range(nt):

        # Apply periodic boundary condition using NumPy roll
        source = [source_term.subs({x: x_pos, t: ti}).evalf() for x_pos in mesh]


        u_new = u0 - coeff * (u0 - np.roll(u0, 1)) + dt * np.array(source)

        u0 = u_new.copy()
        u_all.append(u0)
        ti+=dt;

    return u_all


# Simulation parameters
L = 1.0        # Length of the domain
nx = 201       # Number of spatial points
# Definition of the time step based on CFL
CFL=0.01
c = 1.0        # Advection speed
tfinal=1e-3

# Spatial grid
dx = L / (nx - 1)
mesh = np.linspace(0, L, nx)
dt=CFL*dx/(c)
print("Time step=",dt)
nt=int(tfinal/dt)

x, t = symbols('x t')
u_mms = sin(2*np.pi*x) * cos(np.pi*t)
print("Manufactured solution :",u_mms)

expr=compute_manufactured_solution(u_mms,c, x, t)
u_0=np.sin(2*np.pi*mesh)
u=advect_vectorized(u_0,mesh,dt,nt,c,expr)

tfinal=nt*dt

u_ref=[u_mms.subs({x: x_pos, t: tfinal}).evalf() for x_pos in mesh]

plt.figure(figsize=(10, 6))
plt.plot(mesh, np.array(u)[-1,:],label="solution at "+str(tfinal),color='r')
plt.plot(mesh[::10], np.array(u_ref)[::10],label="u_ref",marker='o',linestyle='None')
plt.legend()
plt.grid(True)
plt.show()
