import mutationpp as mpp
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines   #pip install matplotlib-label-lines

#Exercise 1 - Equilibrium

# Define mixture name",
mixtureName = "air_5_50mbar"
#mixtureName = "air_11"
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
Tout = 8000.0 
#Tout = 20000.0 #for air11
#Tout = 15000 #for CO2_8
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
    plt.plot(
        Temperatures,
        species_descriptor[:, i],
        label=mix.speciesName(i),
        linewidth=2.5  # Épaisseur du trait
    )
    ax = plt.gca()

# Labels sur les lignes
xvals = [6000, 5000, 3000, 3900, 1800]
labelLines(ax.get_lines(), align=True, fontsize=12, ha='left', xvals=xvals)

# Labels d’axes (grande taille)
plt.ylabel("mole fraction", rotation=90, fontsize=16)
plt.xlabel('Temperature, K', fontsize=16)

# Style des ticks
ax.tick_params(
    axis='both',
    labelsize=14,     # taille du texte
    width=2,          # épaisseur des petits traits (ticks)
    length=8,         # longueur des ticks
    direction='inout' # ticks vers l'intérieur ET l'extérieur
)

# Optionnel : rend les ticks visibles de tous côtés
ax.tick_params(top=True, right=True)

# Nettoyage esthétique
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.show()
    
#print("At T = ", Tout, " K")
#for iSp in range(mix.nSpecies()):
    #print("Molar Fraction of ", mix.speciesName(iSp), "is ", mix.X()[iSp])
