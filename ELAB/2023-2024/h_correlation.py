import numpy as np

class const:
    # general constants
    g = 9.81

    # constants relative to experiment
    D = 0.46
    T_inf = 22.5 + 273.15
    d_nozzle = 1.26e-2

    # tabulated values for air properties
    T_prop = (300, 350)
    viscosity = (15.89e-6, 20.92e-6)
    alpha = (22.5e-6, 29.9e-6)
    prandtl = (0.707, 0.700)
    conductivity = (26.3e-3, 30.0e-3)
    density = (1.1614, 0.9950)


def lin_interp(T12, y12, T):
    slope = (y12[1] - y12[0]) / (T12[1] - T12[0])
    return slope*(T - T12[0]) + y12[0]


def h_natural(Tw:float) -> float:
    # air properties at Tf
    Tf = (Tw + const.T_inf)/2
    beta = 1/Tf
    visc = lin_interp(const.T_prop, const.viscosity, Tf)
    alpha = lin_interp(const.T_prop, const.alpha, Tf)
    Pr = lin_interp(const.T_prop, const.prandtl, Tf)
    k = lin_interp(const.T_prop, const.conductivity, Tf)

    # Rayleigh number
    Ra = const.g*beta*(Tw-const.T_inf)*const.D**3 / (visc*alpha)

    # Nusselt number
    Nu = (0.825 + 0.387*Ra**(1/6) / (1 + (0.492/Pr)**(9/16))**(8/27))**2

    # average heat transfer coefficient
    h = k*Nu / const.D

    return h

def h_forced(T2:float, Vdot:float, r:float, dist:float) -> float:
    # air properties at Tf
    Tf = (T2 + const.T_inf)/2
    visc = lin_interp(const.T_prop, const.viscosity, Tf)
    Pr = lin_interp(const.T_prop, const.prandtl, Tf)
    k = lin_interp(const.T_prop, const.conductivity, Tf)

    # correlation parameters
    Ar = const.d_nozzle**2/(4*r**2)
    if Ar < 0.004 or Ar > 0.04:
        print("Careful ! Ar = " + str(Ar) + " out of bounds [0.004, 0.04]")

    H_D = dist/const.d_nozzle
    if H_D < 2 or H_D > 12:
        print("Careful ! H/D = " + str(H_D) + " out of bounds [2, 12]")

    Re = 4*Vdot / (np.pi*visc*const.d_nozzle)
    if Re < 2e3 or Re > 4e5:
        print("Careful ! Re = " + str(Re) + " out of bounds [2e3, 4e5]")

    # correlation
    G = 2*Ar**(0.5) * (1-2.2*Ar**0.5)/(1 + 0.2*(H_D-6)*Ar**0.5)
    Nu = G * (2*Re**0.5*(1 + 0.005*Re**0.55)**0.5) * Pr**0.42
    
    return (k*Nu / const.d_nozzle)

