import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.linalg import expm


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
dt=t_eval[2]-t_eval[1] # delta t.


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


#%% --------------------------------------------------------------------------
# Run the full system as before
#%----------------------------------------------------------------------------


# we run the first 50 seconds:
sol = solve_ivp(turbine_ode, t_span, [1], t_eval=t_eval, method='RK45')

# Extract results
t_ode = sol.t
omega_ode = sol.y[0]
v_w_ode = np.array([wind_speed(ti) for ti in t_ode])
P_g_ode = eta_g * T_g * omega_ode



#%% --------------------------------------------------------------------------
# Kalman Filter estimate
#%----------------------------------------------------------------------------

n_I=200 # index at which we start the assimilation
# Initial state for the assimilation
s0_full_ode=[omega_ode[n_I],v_w_ode[n_I]]
h0_full_ode=P_g_ode[n_I]


n_t=len(t_eval)-n_I # number of samples
t_kalman=t_eval[n_I::]



def f(s):
    """Nonlinear dynamics of the pendulum."""
    omega, v_w = s # 
    ds1 = 1/J*(aerodynamic_torque(v_w, omega)-T_g)
    ds2 = 0
    return np.array([ds1, ds2])

def dCp_dlambda_f(lambda_):
    if lambda_ <= 0:
        raise ValueError("lambda_ must be positive.")
    term1 = 319 / lambda_**3
    term2 = 39.02 / lambda_**2
    exponent = np.exp(-12.5 / lambda_)
    return exponent * (term1 - term2)

def dF_ds(s):
    """Compute the Jacobian of f(s) with respect to the state s."""
    omega, v_w = s
    lambda_=omega*Rad/v_w
    Cp=power_coefficient(lambda_).item() # make sure its just a number with .item()
    dCp_dlambda=dCp_dlambda_f(lambda_).item()
    M=rho*A/(2*J) # multiplicative term
    # first term
    dCp_domega=dCp_dlambda*Rad/v_w
    df_domega=M*v_w**3*(dCp_domega/omega-Cp/omega**2)    
    # second term
    dCp_dv_w=dCp_dlambda*(-omega*Rad/v_w**2)
    df_dv_w=M/omega*(dCp_dv_w*v_w**3+3*Cp*v_w**2)
    
    Jacobian= np.array([[df_domega,df_dv_w],[0,0]])
         
    return Jacobian

def rk4_step(s, dt):
    """Runge-Kutta 4th order integration step."""
    k1 = f(s)
    k2 = f(s + k1 * dt / 2)
    k3 = f(s + k2 * dt / 2)
    k4 = f(s + k3 * dt)
    return s + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6


# =============================================================================
# Initialize State, Covariance and Observation Matrices
# =============================================================================

H = np.array([[eta_g*T_g/1e6, 0]])  # Measure angle only
initial_uncertainty = np.array([[0.01, 0], [0, 1.0]])  # Larger uncertainty in s2
x_0 = np.array([2.9, 10])  # Initialize s1: here you guess!

# =============================================================================
R = np.array([[0.2]])  # Measurement noise for the power (in MW!)
Q = np.array([[0.1, 0], [0, 10]])  # Process noise covariance
# =============================================================================



s_p = np.zeros((2,n_t))  # Predicted state over time
Sigma_p = np.zeros((2, 2, n_t))  # Predicted covariance matrices over time
Sigma_p[:, :, 0] = initial_uncertainty  # Initial uncertainty
s_p[:, 0] = x_0  # Adjusted initial state based on observation mode

delta_s1 = np.zeros(n_t)  # 95% CI for angle
delta_s2 = np.zeros(n_t)  # 95% CI for angular velocity
delta_s1[0] = 1.96 * np.sqrt(Sigma_p[0, 0, 0])
delta_s2[0] = 1.96 * np.sqrt(Sigma_p[1, 1, 0])

#% =============================================================================
# Extended Kalman Filter Propagation
# =============================================================================
for k in range(n_t - 1):
    # Update step --------------------------------------------
    z_k = P_g_ode[n_I+k]/1e6+np.random.normal(loc=0, scale=0.05, size=1)  # Measurement at the current time step
    K = Sigma_p[:, :, k] @ H.T / (H @ Sigma_p[:, :, k] @ H.T + R)  # Scalar case
   
    s_p[:, k] = s_p[:, k] + K @ (z_k - H @ s_p[:, k])  # Correct state estimate
    Sigma_p[:, :, k] = (np.eye(2) - K @ H) @ Sigma_p[:, :, k]  # Correct covariance

    # Prediction step ----------------------------------------
    F_k =dF_ds(s_p[:, k])  # Continuous-time Jacobian
    F_d = expm(F_k * dt)        # Discrete-time Jacobian
    s_p[:, k + 1] = rk4_step(s_p[:, k], dt)  # Nonlinear state prediction
    Sigma_p[:, :, k + 1] = F_d @ Sigma_p[:, :, k] @ F_d.T + Q  # Covariance prediction
    
    # Compute uncertainties (95% CI)
    delta_s1[k + 1] = 1.96 * np.sqrt(Sigma_p[0, 0, k + 1])
    delta_s2[k + 1] = 1.96 * np.sqrt(Sigma_p[1, 1, k + 1])


# Shot the full solution and comparison with truth


# Plot results
plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(t_ode, v_w_ode, label='Wind Speed (m/s)')
plt.plot(t_kalman,s_p[1,:],label='Inferred',color='red')
plt.plot()
plt.ylabel('Wind Speed (m/s)')
plt.legend()
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t_ode, omega_ode, label='Rotor Speed (rad/s)', color='orange')
plt.plot(t_kalman,s_p[0,:],label='Inferred',color='red')

plt.ylabel('Rotor Speed (rad/s)')
plt.legend()
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t_ode, P_g_ode/1e6, label='Generated Power (W)', color='green')
plt.plot(t_kalman,eta_g*T_g*s_p[0,:]/1e6,label='Inferred',color='red')
plt.xlabel('Time (s)')
plt.ylabel('Generated Power (kW)')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()






