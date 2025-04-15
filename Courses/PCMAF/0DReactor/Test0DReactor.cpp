#include "/home/jpe/Mutationpp/src/general/mutation++.h"
#include <Eigen/Dense>

using namespace Mutation;
using namespace Eigen;
using namespace std;

// Fonction compatible C pour appel Python
extern "C" int ZeroDReactor(const char* mixture_file, double Tinlet, double Tsurface, double Pstat) {
    // Initialisation du gaz
    MixtureOptions opts(mixture_file);
    Mixture mix(opts);        

    // General settings
    const size_t set_state_with_rhoi_T = 1;
    const size_t pos_T_trans = 0;

    // Number of species
    size_t ns = mix.nSpecies();
    // Number of energy equations
    size_t nT = mix.nEnergyEqns();

    std::cout << "\n--- Gas ---" << std::endl;

    // Equilibrium computation
    mix.equilibrate(Tinlet, Pstat);

    std::cout << "Species mass fractions:\n";
    for (int i = 0; i < mix.nSpecies(); ++i)
        std::cout << mix.speciesName(i) << ": " << mix.Y()[i] << '\n';     

    // Computation of the species densities
    VectorXd rhoi_s(ns);   
    mix.densities(rhoi_s.data());

    std::cout << "\n--- Surface ---" << std::endl;
    
    // Computation of species mass fraction
    VectorXd xi_e(ns);
    xi_e = Map<const VectorXd>(mix.X(), ns);
    double dx = 1e-3; 


    // Computation of surface mass balance
    mix.setSurfaceState(rhoi_s.data(), &Tsurface, set_state_with_rhoi_T);
    mix.setDiffusionModel(xi_e.data(), dx);
    mix.setIterationsSurfaceBalance(3000);
    mix.solveSurfaceBalance();
    mix.getSurfaceState(rhoi_s.data(), &Tsurface, set_state_with_rhoi_T);

    std::cout << "\nSurface mass fractions: " << std::endl;
    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        std::cout << mix.speciesName(i) << ": " << mix.Y()[i] << std::endl;
    }

    // Use the results of surface mass balance to compute reaction rates
    mix.setState(rhoi_s.data(), &Tsurface, set_state_with_rhoi_T);

    // Computation of reaction rates
    VectorXd wdot(ns);
    mix.surfaceReactionRates(wdot.data());

    std::cout << "\nSurface reaction rates: " << std::endl;
    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        std::cout << mix.speciesName(i) << ": " << wdot.data()[i] << '\n';
    }

    // Computation of thermal conductivity
    double lambda = mix.equilibriumThermalConductivity();
    std::cout << "Lambda = " << lambda << std::endl;

    // Computation of conduction heat flux
    double qw = (Tinlet - Tsurface) / dx * lambda;
    std::cout << "qw = " << qw << std::endl;

    // Computation of unit mass species enthlapy
    VectorXd h_s(ns);
    mix.getEnthalpiesMass(h_s.data());

    // Computation of species enthlapy
    for (int i = 0 ; i < ns; i ++)
    {
        h_s[i] = h_s[i] * mix.Y()[i];
    }

    // Computation of total heat flux (conduction + catalysis)
    for (int i = 0; i < ns; i++) {
        qw += wdot[i] * h_s[i];
    }

    std::cout << "qw = " << qw << std::endl;

    return 0;
}
