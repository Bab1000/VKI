import numpy as np
from scipy.optimize import root_scalar

def LoadNASACoeffs():

    NASAcoefficients = {
        'O': {
            'high_T': np.array([1.779004264E+08, -1.082328257E+05, 2.810778365E+01, -2.975232262E-03, 1.854997534E-07, -5.796231540E-12, 7.191720164E-17, 8.890942630E+05]),

            'mid_T': np.array([2.619020262E+05, -7.298722030E+02, 3.317177270E+00, -4.281334360E-04, 1.036104594E-07, -9.438304330E-12, 2.725038297E-16, 3.392428060E+04]),

            'low_T':  np.array([-7.953611300E+03, 1.607177787E+02, 1.966226438E+00, 1.013670310E-03, -1.110415423E-06, 6.517507500E-10, -1.584779251E-13, 2.840362437E+04])

        },
        'O2': {
            'high_T': np.array([4.975294300E+08, -2.866106874E+05, 6.690352250E+01, -6.169959020E-03, 3.016396027E-07, -7.421416600E-12, 7.278175770E-17, 2.293554027E+06]),

            'mid_T': np.array([-1.037939022E+06, 2.344830282E+03, 1.819732036E+00, 1.267847582E-03, -2.188067988E-07, 2.053719572E-11, -8.193467050E-16, -1.689010929E+04]),

            'low_T':  np.array([-3.425563420E+04, 4.847000970E+02, 1.119010961E+00, 4.293889240E-03, -6.836300520E-07, -2.023372700E-09, 1.039040018E-12, -3.391454870E+03])

        },
        'N': {
            'high_T': np.array([5.475181050E+08, -3.107574980E+05, 6.916782740E+01, -6.847988130E-03, 3.827572400E-07, -1.098367709E-11, 1.277986024E-16, 2.550585618E+06]),

            'mid_T': np.array([8.876501380E+04, -1.071231500E+02, 2.362188287E+00, 2.916720081E-04, -1.729515100E-07, 4.012657880E-11, -2.677227571E-15, 5.697351330E+04]),

            'low_T':  np.array([0.000000000E+00, 0.000000000E+00, 2.500000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 5.610463780E+04])

        },
        'NO': {
            'high_T': np.array([-9.575303540E+08, 5.912434480E+05, -1.384566826E+02, 1.694339403E-02, -1.007351096E-06, 2.912584076E-11, -3.295109350E-16, -4.677501240E+06]),

            'mid_T': np.array([2.239018716E+05, -1.289651623E+03, 5.433936030E+00, -3.656034900E-04, 9.880966450E-08, -1.416076856E-11, 9.380184620E-16, 1.750317656E+04]),

            'low_T':  np.array([-1.143916503E+04, 1.536467592E+02, 3.431468730E+00, -2.668592368E-03, 8.481399120E-06, -7.685111050E-09, 2.386797655E-12, 9.098214410E+03])

        },
        'N2': {
            'high_T': np.array([8.310139160E+08, -6.420733540E+05, 2.020264635E+02, -3.065092046E-02, 2.486903333E-06, -9.705954110E-11, 1.437538881E-15, 4.938707040E+06]),

            'mid_T': np.array([5.877124060E+05, -2.239249073E+03, 6.066949220E+00, -6.139685500E-04, 1.491806679E-07, -1.923105485E-11, 1.061954386E-15, 1.283210415E+04]),

            'low_T':  np.array([2.210371497E+04, -3.818461820E+02, 6.082738360E+00, -8.530914410E-03, 1.384646189E-05, -9.625793620E-09, 2.519705809E-12, 7.108460860E+02])

        }
    }

    return NASAcoefficients

def LoadMolarMass():

    SpeciesMolarMasses = {

        'O': 15.9994e-3, # kg/mol
        'O2': 32e-3,
        'N': 14.0067e-3, 
        'NO': 30.01e-3,
        'N2': 28.0134e-3
    }

    return SpeciesMolarMasses

def NASA_cp_species(species,T):

    # Loading the coefficients
    NASAcoefficients = LoadNASACoeffs()

    # Getting the coefficients
    if T >= 6000:
        coeffs = NASAcoefficients[species]['high_T']
    elif T >= 1000 and T < 6000:
        coeffs = NASAcoefficients[species]['mid_T']
    elif T < 1000:
        coeffs = NASAcoefficients[species]['low_T']

    # Computing Cpi/Ri
    cpi_Ri = coeffs[0]*T**-2 + coeffs[1]*T**-1 + coeffs[2] + coeffs[3]*T + coeffs[4]*T**2 + coeffs[5]*T**3 + coeffs[6]*T**4

    # Loading species molar mass
    SpeciesMolarMass = LoadMolarMass()

    # Gathering species molar mass
    MM_i = SpeciesMolarMass[species]
    
    # Computing Ri 
    R_univ = 8.314462618 # J/mol/K
    Ri = R_univ / MM_i  # J/kg/K

    # Computing species Cp
    cpi = cpi_Ri * Ri

    return cpi,Ri

def gamma_T(cp,R):
    # Computing cv
    cv = cp - R
    # Computing species gamma
    gamma_i = cp / cv
    return gamma_i

def LoadMassFraction():
    MassFraction = {
        "N" : 2.39879e-80,
        "O" : 1.25697e-41,
        "NO" : 2.69731e-16,
        "N2" : 0.767082,
        "O2" : 0.232918
    }
    return MassFraction

def stagnation_enthalpy_residual(M1,H0_target,list_species,T,gamma):

    # Storage for species enthalpy
    h_i_list = []
    R_i_list = []
    Cp_i_list = []

    # Looping over all the species in the air mixture
    for species in list_species:

        # Load coefficients
        NASAcoefficients = LoadNASACoeffs()

        # Getting the coefficients
        if T >= 6000:
            coeffs = NASAcoefficients[species]['high_T']
        elif T >= 1000 and T < 6000:
            coeffs = NASAcoefficients[species]['mid_T']
        elif T < 1000:
            coeffs = NASAcoefficients[species]['low_T']

        #Computing cp
        cp_i,Ri = NASA_cp_species(species,T)
        # Computing species enthalpy
        h_i = cp_i * T
        #h_i_RiT = -coeffs[0]*T**-2 + coeffs[1]*np.log(T)*T**-1 + coeffs[2] + coeffs[3]*T/2 + coeffs[4]*T**2/3 + coeffs[5]*T**3/4 + coeffs[6]*T**4/5 + coeffs[6]/T
        #h_i = h_i_RiT * Ri * T
        # Saving species enthalpy
        h_i_list.append(h_i)
        # Saving species cp
        Cp_i_list.append(cp_i)
        # Saving species R
        R_i_list.append(Ri)

    # Loading species molar mass
    MassFraction = LoadMassFraction()
    MolarMass = LoadMolarMass()

    R = 0.0
    cp = 0.0
    h = 0.0

    for i in range(len(list_species)):
        # Gathering the species
        species = list_species[i]
        # Gathering the species molar mass
        Y_i = MassFraction[species]
        #Gathering species Molar mass
        MM_i = MolarMass[species]

        # Computing global R
        R += Y_i * R_i_list[i]
        # Computing global cp
        cp += Y_i * Cp_i_list[i]
        # Computing global h
        h += Y_i * h_i_list[i]

    # Computing gamma
    #gamma = gamma_T(cp,R)

    V_squared = M1**2 * gamma * R * T

    return h + 0.5 * V_squared - H0_target

def UpstreamT(M1,H0_target,list_species,gamma):
    # Solve for T1
    sol = root_scalar(
        lambda T: stagnation_enthalpy_residual(M1,H0_target,list_species,T,gamma),
        bracket=[200, 20000],
        method='brentq'
    )
    
    T1_solution = sol.root if sol.converged else None
    
    print(f"Static Temperature upstream of shock (air5): {sol.root:.2f} K")   

    return T1_solution 