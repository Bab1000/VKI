import mutationpp as mpp
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines   #pip install matplotlib-label-lines

#Exercise 1 - Equilibrium

# Define mixture name",
#mixtureName = "air_5"
mixtureName = "air_11_50mbar"
#mixtureName = "CO2_8"
myMixtureOptions = mpp.MixtureOptions(mixtureName)
print("Mixture species:", myMixtureOptions.getSpeciesDescriptor())

#Define mixture State Model
stateModelName = "Equil"
myMixtureOptions.setStateModel(stateModelName)
print("Mixture state model:", myMixtureOptions.getStateModel())

#Create Mixture
mix = mpp.Mixture(myMixtureOptions)

#Compute Equilibrium composition at p = 1 atm and T in [500, 8000] K 
ns = mix.nSpecies()
Tin = 500.0
#Tout = 8000.0 
Tout = 15000 #for CO2_8 and Air11
Temperatures = np.linspace(Tin, Tout, 100)
Pressure = 101325

#mix.addComposition("N2:0.79, O2:0.21", True) #default
#mix.addComposition("N2:0.5, O2:0.5", True)

species_descriptor = np.empty((0, ns), float)
for Temperature in Temperatures:
	mix.setState(Pressure, Temperature, 1)
	species_descriptor = np.vstack((species_descriptor, np.array(mix.X())))
	if Temperature == Tin:
		print("At T = ", Temperature, " K")
		for iSp in range(mix.nSpecies()):
			print("Molar Fraction of ", mix.speciesName(iSp), "is ", mix.X()[iSp])
		print("\n")

plt.figure(figsize=(12,8))
plt.tight_layout()
for i in range(ns):
    plt.plot(Temperatures, species_descriptor[:, i], label=mix.speciesName(i))
    ax = plt.gca()

    #xvals = [6000, 5000, 3000, 3900, 1800] #position of the label Air 5
    xvals = [15000, 14000, 13000, 12000, 16000, 17000, 6000, 5000, 3000, 3900, 1800] #position of the label Air 11
    #xvals = [13500, 13000, 14000, 7900, 7000, 5000, 700, 3200]

    labelLines(plt.gca().get_lines(), align=True, fontsize=10, ha='left', xvals=xvals)

    plt.ylabel("mole fraction", rotation=90)
    plt.xlabel('Temperature, K')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

plt.show()
    
#print("At T = ", Tout, " K")
#for iSp in range(mix.nSpecies()):
    #print("Molar Fraction of ", mix.speciesName(iSp), "is ", mix.X()[iSp])


iT = 0
viscosity = np.zeros(100, float)
thermalConductivity = np.zeros((100, 4), float)

for Temperature in Temperatures:
    mix.setState(Pressure, Temperature, 1)
    
    viscosity[iT] = mix.viscosity()
    thermalConductivity[iT, 0] = mix.heavyThermalConductivity()
    thermalConductivity[iT, 1] = mix.internalThermalConductivity(Temperature)
    thermalConductivity[iT, 2] = mix.reactiveThermalConductivity()
    thermalConductivity[iT, 3] = mix.heavyThermalConductivity() + mix.internalThermalConductivity(Temperature)+ mix.reactiveThermalConductivity()
    iT+=1

#viscosity    
plt.figure()
plt.plot(Temperatures, viscosity, color='k')
ax = plt.gca()
plt.ylabel(r'$\eta$ [Pa s]', rotation=90)
plt.xlabel('Temperature, K')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()

#plt.savefig(\"viscosityair11.png\", format='png', transparent=True)\n",

plt.figure()
listColor = ["b", "g", "r", "k"]
listLabel = ["Heavy", "Internal", "Reactive", "Total"]
for j in range(4):
    plt.plot(Temperatures, thermalConductivity[:, j], color = listColor[j] )
ax = plt.gca()
plt.ylabel(r'$\lambda$ [W/mK]', rotation=90)
plt.xlabel('Temperature, K')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.text(6000, 3.3, "Heavy", color = "k")
plt.text(6700, 2.5, "Reactive", color = "r")
plt.text(6700, 0.5, "Heavy", color = "b")
plt.text(6800, .18, "Internal", color = "g")
plt.show()
    
