# Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import pandas as pd
from scipy.integrate import odeint


# =============================================================================
# Load and Preprocess the Measurement Data
# =============================================================================

# Measurement file number
MEAS = 1
file_name = f'data_pendulum/measurement_{MEAS}.dat'

# Load measurement data
data = pd.read_csv(file_name, delimiter='\t', skiprows=1)
t_s = data.iloc[:, 0].values  # Time vector
s1_s = data.iloc[:, 1].values  # Measured angle (theta)
s2_s = data.iloc[:, 2].values  # Measured angular velocity (theta_dot)
s_s = np.array([s1_s, s2_s])  # Combined state vector

# Extract parameters
n_t = len(t_s)  # Number of time steps
dt = t_s[1] - t_s[0]  # Time step

# =============================================================================
# Define System Parameters and Functions
# =============================================================================
# we compute also the true solution:
m=120/1000 # mass of the pendulum  [kg]
L=41/100 # full length of the pendulum  [m]
w=37/1000 # width of the pendulum [m]
l_cm=24/100 # distance from center of mass to pivot [m]
g=9.815 # gravitational acceleration [m/s^2]
I_x=1/12*m*(L**2+w**2)+m*l_cm**2 # Moment of Inertia [kg m^2]

omega_n=np.sqrt(m*g*l_cm/I_x) # natural omega
degs=np.array([+35,-28,+15,-45,36,-22,18,42])
s0_s=degs*np.pi/180
mu_s=np.array([0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18])

mu=mu_s[MEAS]      
s0=[s0_s[MEAS],0]


def f(s):
    """Nonlinear dynamics of the pendulum."""
    s1, s2 = s
    ds1 = s2
    ds2 = -mu * s2 - omega_n**2 * np.sin(s1)
    return np.array([ds1, ds2])

def compute_F(s):
    """Compute the Jacobian of f(s) with respect to the state s."""
    s1, _ = s
    return np.array([
        [0, 1],
        [-omega_n**2 * np.cos(s1), -mu]
    ])

def rk4_step(s, dt):
    """Runge-Kutta 4th order integration step."""
    k1 = f(s)
    k2 = f(s + k1 * dt / 2)
    k3 = f(s + k2 * dt / 2)
    k4 = f(s + k3 * dt)
    return s + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6

#%% =============================================================================
# Define Observation Options
# =============================================================================
# Select which state(s) to measure: "s1", "s2", or "both"
observed_state = "both"  # Change to "s1", "s2", or "both"

if observed_state == "s1":
    H = np.array([[1, 0]])  # Measure angle only
    R = np.array([[0.002]])  # Measurement noise for angle
    initial_uncertainty = np.array([[0.001, 0], [0, 1.0]])  # Larger uncertainty in s2
    #x_0 = np.array([s1_s[0], 0.0])  # Initialize s1 from measurement, s2 assumed zero
    # or take the true ones (for MEAS=1)
    x_0=s0
    Q = np.array([[0.002, 0], [0, 0.1]])  # Process noise covariance
    z_k_func = lambda k: np.array([s1_s[k]])  # Measurement vector
elif observed_state == "s2":
    H = np.array([[0, 1]])  # Measure angular velocity only
    R = np.array([[0.08]])  # Measurement noise for angular velocity
    initial_uncertainty = np.array([[0.01, 0], [0, 0.1]])   # Small initial uncertainty
    Q = np.array([[0.002, 0], [0, 0.08]])  # Process noise covariance
    #x_0 = np.array([-0.45, s2_s[0]])  # Initialize s2 from measurement, s1 assumed -0.45
    # or take the true ones (for MEAS=1)
    x_0=s0 
    z_k_func = lambda k: np.array([s2_s[k]])  # Measurement vector
elif observed_state == "both":
    H = np.eye(2)  # Measure both states
    R = np.array([[0.002, 0], [0, 0.08]])  # Measurement noise for both states
    initial_uncertainty = np.array([[0.01, 0], [0, 0.1]])   # Small initial uncertainty
    Q = R/10000  # Process noise covariance (trust more the model or the measurements?)
    x_0 = np.array([s1_s[0], s2_s[0]])  # Initialize both states from measurements
    # or take the true ones (for MEAS=1)
    #x_0=s0    
    z_k_func = lambda k: s_s[:, k]  # Measurement vector
else:
    raise ValueError("Invalid observed_state. Choose 's1', 's2', or 'both'.")

# =============================================================================
# Initialize State and Covariance Matrices
# =============================================================================

s_p = np.zeros_like(s_s)  # Predicted state over time
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
    z_k = z_k_func(k)  # Measurement at the current time step
    if R.size == 1:
        K = Sigma_p[:, :, k] @ H.T / (H @ Sigma_p[:, :, k] @ H.T + R)  # Scalar case
    else:
        K = Sigma_p[:, :, k] @ H.T @ np.linalg.inv(H @ Sigma_p[:, :, k] @ H.T + R)  # Matrix case

    s_p[:, k] = s_p[:, k] + K @ (z_k - H @ s_p[:, k])  # Correct state estimate
    Sigma_p[:, :, k] = (np.eye(2) - K @ H) @ Sigma_p[:, :, k]  # Correct covariance

    # Prediction step ----------------------------------------
    F_k = compute_F(s_p[:, k])  # Continuous-time Jacobian
    F_d = expm(F_k * dt)        # Discrete-time Jacobian
    s_p[:, k + 1] = rk4_step(s_p[:, k], dt)  # Nonlinear state prediction
    Sigma_p[:, :, k + 1] = F_d @ Sigma_p[:, :, k] @ F_d.T + Q  # Covariance prediction
    # Compute uncertainties (95% CI)
    delta_s1[k + 1] = 1.96 * np.sqrt(Sigma_p[0, 0, k + 1])
    delta_s2[k + 1] = 1.96 * np.sqrt(Sigma_p[1, 1, k + 1])

# =============================================================================
# Visualization of Results
# =============================================================================


    
# We define the function of the system    
def f_ODEINT(s, t,mu=0.5,omega_n=5):
    # this is the function f in the slides
    s_1_dot = s[1]
    s_2_dot = -mu*s[1]-omega_n**2*np.sin(s[0])
    return [s_1_dot, s_2_dot]


Y = odeint(f_ODEINT, s0, t_s,args=(mu,omega_n))
s1_f=Y[:,0]; s2_f=Y[:,1] # Get the solution





# Create a single plot with subfigures for s1 (angle) and s2 (angular velocity)
fig, axs = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

# Title for the entire figure based on the observed_state
fig.suptitle(f"EKF Predictions: Observing {observed_state.upper()}", fontsize=18)

# Subfigure 1: Plot for angle (s1)
axs[0].plot(t_s, s_p[0, :], label="Predicted $s_1$ (Angle)")
axs[0].plot(t_s, s1_s, 'r--', label="Measured $s_1$")
axs[0].plot(t_s,s1_f,'k',label="Truth")

axs[0].fill_between(t_s,
                    s_p[0, :] - delta_s1,
                    s_p[0, :] + delta_s1,
                    color='blue', alpha=0.2, label="95% CI for $s_1$")
axs[0].set_xlabel("$t$")
axs[0].set_ylabel("$\\theta(t)$ (Angle)")
axs[0].set_title("Angle Prediction ($s_1$)", fontsize=16)
axs[0].grid()
axs[0].legend()

# Subfigure 2: Plot for angular velocity (s2)
axs[1].plot(t_s, s_p[1, :], label="Predicted $s_2$ (Angular Velocity)")
axs[1].plot(t_s, s2_s, 'r--', label="Measured $s_2$")
axs[1].plot(t_s,s2_f,'k',label="Truth")

axs[1].fill_between(t_s,
                    s_p[1, :] - delta_s2,
                    s_p[1, :] + delta_s2,
                    color='blue', alpha=0.2, label="95% CI for $s_2$")
axs[1].set_xlabel("$t$")
axs[1].set_ylabel("$\dot{\\theta}(t)$ (Angular Velocity)")
axs[1].set_title("Angular Velocity Prediction ($s_2$)", fontsize=16)
axs[1].grid()
axs[1].legend()

# Show the combined plot
plt.show()

output_filename = f"EKF_Predictions_Observing_{observed_state.upper()}.png"

fig.savefig(output_filename, dpi=300)
