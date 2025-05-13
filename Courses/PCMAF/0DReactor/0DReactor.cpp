#include "mutation++.h"
#include <Eigen/Dense>

using namespace Mutation;
using namespace Eigen;
using namespace std;

// Fonction compatible C pour appel Python
extern "C" int ZeroDReactor(const char* mixture_file, double Tinlet, double Tsurface, double Pstat, double dx, double* wdot_out, double* qw_out ) {
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

    // Equilibrium computation
    mix.equilibrate(Tinlet, Pstat);

    // Computation of the species densities
    VectorXd rhoi_s(ns);   
    mix.densities(rhoi_s.data());
    
    // Computation of species mass fraction
    VectorXd xi_e(ns);
    xi_e = Map<const VectorXd>(mix.X(), ns);

    // Computation of surface mass balance
    double dx_diff = 1e-3;
    mix.setSurfaceState(rhoi_s.data(), &Tsurface, set_state_with_rhoi_T);
    mix.setDiffusionModel(xi_e.data(), dx_diff);
    mix.setIterationsSurfaceBalance(3000);
    mix.solveSurfaceBalance();
    mix.getSurfaceState(rhoi_s.data(), &Tsurface, set_state_with_rhoi_T);

    // Use the results of surface mass balance to compute reaction rates
    mix.setState(rhoi_s.data(), &Tsurface, set_state_with_rhoi_T);

    // Computation of reaction rates
    VectorXd wdot(ns);
    mix.surfaceReactionRates(wdot.data());
    //std::cout << "\nSurface reaction rates: " << std::endl;
    //for (int i = 0; i < mix.nSpecies(); ++i)
    //{
    //    std::cout << mix.speciesName(i) << ": " << wdot.data()[i] << '\n';
    //}

    // Computation of thermal conductivity
    double lambda = mix.equilibriumThermalConductivity();

    // Computation of conduction heat flux
    //std::cout<<"---> [INFO] Thermal Conductivity : "<< lambda << endl;
    double qw = - ((Tinlet - Tsurface) / dx) * lambda;

    //std::cout<<"---> [INFO] Conduction heat flux : "<< qw << endl;

    // Computation of unit mass species enthlapy
    VectorXd h_s(ns);

    mix.speciesHOverRT(Tinlet,h_s.data());

    // Computation of species enthlapy
    for (int i = 0 ; i < ns; i++)
    {
        h_s[i] = h_s[i] * 8.314411 * Tinlet;
    }

    // Computation of total heat flux (conduction + catalysis)
    for (int i = 0; i < ns; i++) 
    {
        qw += wdot[i] * h_s[i];
    }

    // Cpying the values to return it in the python world
    *qw_out = std::abs(qw);
    for (int i = 0; i < ns; i++) 
    {
        wdot_out[i] = wdot[i];  
    }

    //std::cout<<"qw_out : "<< std::abs(qw) <<endl;
    //std::cout<<"dx : "<< dx <<endl;

    return 0;
}
