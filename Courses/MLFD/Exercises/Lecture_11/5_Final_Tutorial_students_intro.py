import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Customize plot settings for LaTeX and larger fonts
plt.rcParams.update({
    'text.usetex': True,
    'font.size': 18,
    'font.family': 'serif'
})


#%% Constants
rho = 1.2  # Air density (kg/m^3)
Rad = 40       # Rotor radius (m)
A = np.pi * Rad**2  # Rotor swept area (m^2)
J = 4000000   # Rotor inertia (kg·m^2)
eta_g = 0.98 # Generator efficiency
T_g = 400000   # Generator torque (N·m)

# Time array
t_span = (0, 400)  # Time span for the simulation (seconds)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# Seed for reproducibility
np.random.seed(42)

# Wind speed function
def wind_speed(t):
    term1=10
    term2=3 * np.sin(0.1 * t)*np.exp(-0.001*(t-160)**2)
    term3=2 * np.sin(0.4 * t)*np.exp(-0.005*(t-250)**2)
    term4=0.8 * np.random.normal()
    v_w=term1+term2+term3+term4
    return v_w

# Power coefficient function
def power_coefficient(lambda_):
    if lambda_ <= 0:
        return 0
    return 0.22 * (116 / lambda_ - 5) * np.exp(-12.5 / lambda_)

# Plot the power coefficient:
    
lambdas=np.linspace(2,20,500)    
C_p_F=np.vectorize(power_coefficient)

fig, ax = plt.subplots(figsize=(6,4)) # Create Signal Noisy 
plt.plot(lambdas,C_p_F(lambdas),'ko')
plt.xlabel('$\lambda=\omega R/ v_w$',fontsize=18)
plt.ylabel('$C_p(\lambda)$',fontsize=18)
plt.tight_layout()
plt.ylim(0,0.55)
plt.show()
Name='C_P_Curve.png'
plt.savefig(Name, dpi=300)      
# plt.close(fig)
    

# Aerodynamic torque function
def aerodynamic_torque(v_w, omega):
    if omega == 0:
        return 0
    lambda_ = omega * Rad / v_w
    Cp = power_coefficient(lambda_)
    Pa = 0.5 * rho * A * Cp * v_w**3
    return Pa / omega

# ODE system
def turbine_ode(t, omega):
    v_w = wind_speed(t)
    T_a = aerodynamic_torque(v_w, omega)
    domega_dt = (T_a - T_g) / J
    return domega_dt

# Initial rotor speed
omega_0 = [1.0]  # rad/s
# Solve ODE
sol = solve_ivp(turbine_ode, t_span, omega_0, t_eval=t_eval, method='RK45')

# Calculate power output
omega_sol = sol.y[0]
v_w_sol = np.array([wind_speed(t) for t in sol.t])
P_g = eta_g * T_g * omega_sol

# Plot results
plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(sol.t, v_w_sol, label='Wind Speed (m/s)')
plt.ylabel('Wind Speed (m/s)')
plt.legend()
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(sol.t, omega_sol, label='Rotor Speed (rad/s)', color='orange')
plt.ylabel('Rotor Speed (rad/s)')
plt.legend()
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(sol.t, P_g/1000, label='Generated Power (W)', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Generated Power (kW)')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()



# check the performances-- is that a good torque?
plt.figure(12)
plt.plot(sol.t,omega_sol*Rad/v_w_sol,'ko')




