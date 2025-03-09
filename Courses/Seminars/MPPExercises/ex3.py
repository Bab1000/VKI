import mutationpp as mpp
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines   #pip install matplotlib-label-lines

#Define mixture name
myMixtureOptions = mpp.MixtureOptions() # Empty options object 
myMixtureOptions.setSpeciesDescriptor("O O2") # Setting species
print("Mixture species:", myMixtureOptions.getSpeciesDescriptor())

#Define mixture State Model 
stateModelName = "ChemNonEq1T"
myMixtureOptions.setStateModel(stateModelName) #Equilibrium model
myMixtureOptions.setThermodynamicDatabase("RRHO")
print("Mixture state model:", myMixtureOptions.getStateModel())

#Define Chemical Mechanism
myMixtureOptions.setMechanism("O2")
print("Mixture Chemical Mechanism:", myMixtureOptions.getMechanism())
#Create Mixture
mix = mpp.Mixture(myMixtureOptions)

nS   = mix.nSpecies()
rhoi = np.zeros(nS) # Densities
wdot = np.zeros(nS) # Chemical rate of production

# Set Initial conditions
Tinit    = 10000 # [K]
rho_init = 3.9e-2 # [kg/m^3] (Chosen to have P_init = 1 atm)
rhoi[mix.speciesIndex("O")] = 0.0
rhoi[mix.speciesIndex("O2")] = rho_init
mix.setState(rhoi, Tinit, 1)
print("Mixture Pressure: ", mix.P())

y_vect = mix.Y()
thermalConductivity = mix.frozenThermalConductivity()
T = Tinit

etot = mix.mixtureEnergyMass()*mix.density()
print("Total Energy (Conserved thorugh simulation) = ", etot, " J/m^3")

time = 0.0
tt = time
tout = 1e-9
tfinal = 1e-4
tsampling = 100
k = 0

# Simulation
while time < tfinal:
    # Get the species production rates
    wdot = mix.netProductionRates()
    # Compute \"stable\" timestep based on maximum allowed change in species densities 
    dt =0.0 
    for i in range(nS): 
        dt += 5e-6 *np.abs(rho_init*mix.Y()[i]/wdot[i]) 
    dt = min(dt / nS, tfinal - time) 

    # Integrate in time 
    for i in range(nS): 
        rhoi[i] += wdot[i]*dt 

    time += dt 
    # Update state 
    mix.setState(rhoi, etot, 0) 
    if k % tsampling == 0 or np.abs(time-tfinal)<dt: 
        y_vect = np.vstack((y_vect, mix.Y())) 
        thermalConductivity = np.append(thermalConductivity, mix.frozenThermalConductivity()) 
        T = np.append(T, mix.T()) 
        tt = np.append(tt, time) 
        tout += tout 
        # Display advancement 
        a = round((np.abs(time))*100 / tfinal, 2) 
        print('\r', r"Computing {:.2f}%".format(a), end='')     
    k+=1 

plt.figure(figsize=(12,8)) 
plt.semilogx(tt, y_vect[:, 0], color = "k") 
plt.semilogx(tt, y_vect[:, 1], color = "r") 
ax = plt.gca() 
plt.ylabel(r'Mass Fraction', rotation=90) 
plt.xlabel('Time, s') 
ax.spines['right'].set_visible(False) 
ax.spines['top'].set_visible(False) 
     
plt.text(1e-5, 0.65, "O2", color = "r") 
plt.text(1e-5, 0.35, "O", color = "k") 
     
plt.figure(figsize=(12,8)) 
plt.semilogx(tt, thermalConductivity) 
ax = plt.gca() 
plt.ylabel(r'$\lambda$, [W/mK]', rotation=90) 
plt.xlabel('Time, s') 
ax.spines['right'].set_visible(False) 
ax.spines['top'].set_visible(False) 
plt.show()
     
plt.figure(figsize=(12,8)) 
plt.semilogx(tt, T) 
ax = plt.gca() 
plt.ylabel(r'Temperature, K', rotation=90) 
plt.xlabel('Time, s') 
ax.spines['right'].set_visible(False) 
ax.spines['top'].set_visible(False)
plt.show()
