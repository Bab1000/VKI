#include "StagLineSol.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using std::setw;

/**
 * Generates the location specific data for input in the radiation solver.
 */
class RadiationProcessor
{
public:

    // Create the mesh based on maximum change in temperature between points
    static void process(StagLineSol& s)
    {
        const double min_dist = 0.00175;
        const double max_dist = 0.00805;
        const double max_change_T = 100.0;
    
        std::ofstream f("rad_input.dat");
        
        int width = 20;
        f.precision(10);
        f << std::scientific;
        
        // Interpolate to get wall temperatures


        // Write the header
        f << "Geometry type\n";
       // f << "SPHERICAL\n";
        f << "CARTESIAN\n";
        f << "Boundary conditions: eps1, Tw1(K), eps2, Tw2(K)\n";
        f << "0.85 " << s.getMixAtFacePoint(0).T() << " ";
        f << "1.0 "  << s.getMixAtFacePoint(s.nFaces()-1).T() << "\n";
        f << setw(width) << "#Loc";
        f << setw(width) << "Tr";
        f << setw(width) << "Tv";
        f << setw(width) << "P";
        
        Mutation::Mixture& mix = s.getMixAtCellPoint(0);
        for (int i = 0; i < mix.nSpecies(); ++i)
            f << setw(width) << mix.speciesName(i);
        f << "\n";
        
        int i = 0;
        while (i < s.nVolumes()) {
            // Write the data for this point
            Mutation::Mixture& mix = s.getMixAtCellPoint(i);
            f << setw(width) << s.faceLoc(i) * 100.0;
            f << setw(width) << mix.T();
            f << setw(width) << mix.Tv();
            f << setw(width) << mix.P();
            
            for (int k = 0; k < mix.nSpecies(); ++k)
                f << setw(width) << mix.X()[k];
            f << "\n";
            
            // Find the next point
            int k = i;
            // Satisfy minimum distance between points
            //while (s.faceLoc(k) - s.faceLoc(i) < min_dist && k < s.nVolumes()) 
           //     k++;
            // Satisfy minimum change and maximum distance between points
            double Ti = mix.T();
            while (s.faceLoc(k) - s.faceLoc(i) < max_dist && k < s.nVolumes()) {
                Mutation::Mixture& mix = s.getMixAtCellPoint(k);
                if (std::abs(Ti - mix.T()) > max_change_T) break;
                k++;
            }
            
            i = k;
        }
        
        // Print the last point
        f << setw(width) << s.faceLoc(s.nVolumes()) * 100.0 << "\n";
        f.close();
    }

    /*
    static void process(StagLineSol& s) 
    {
        std::ofstream f("rad_input.dat");
        
        // Find the thermal boundary layer edge
        int ibl;
        for (ibl = 0; ibl < s.nVolumes(); ++ibl)
            if (s.getMixAtCellPoint(ibl+1).T()-s.getMixAtCellPoint(ibl).T() < 5.0) break;
        ibl += 2;
       // ibl = 100;
        cout << ibl << endl;

        // Find the shock non equilibrium boundary
        int ine;
        for (ine = ibl+1; ine < s.nVolumes(); ++ine)
            if (s.getMixAtCellPoint(ine+1).T()-s.getMixAtCellPoint(ine).T() > 5.0) break;
        ine -= 5;
        cout << ine << endl;

        // Find the freestream
        int ifs;
        for (ifs = ine+6; ifs < s.nVolumes()-1; ++ifs)
            if (s.getMixAtCellPoint(ifs+1).T()-s.getMixAtCellPoint(ifs).T() < 1.0) break;
        ifs += 10;
        cout << ifs << endl;

        int width = 20;
        f.precision(10);
        f << std::scientific;
        
std::cout<<"INSIDE radiaion"<<std::endl;

        // Write the header
        f << setw(width) << "#Loc";
        f << setw(width) << "Tr";
        f << setw(width) << "Tv";
        f << setw(width) << "P";
        
        Mutation::Mixture& mix = s.getMixAtCellPoint(0);
        for (int i = 0; i < mix.nSpecies(); ++i)
            f << setw(width) << mix.speciesName(i);
        f << "\n";
        
        // Now write the location data
        int i = 0;
        while (i < s.nVolumes()) {
            Mutation::Mixture& mix = s.getMixAtCellPoint(i);
            f << setw(width) << s.faceLoc(i) * 100.0;
            f << setw(width) << mix.T();
            f << setw(width) << mix.Tv();
            f << setw(width) << mix.P();
            
            for (int k = 0; k < mix.nSpecies(); ++k)
                f << setw(width) << mix.X()[k];
            f << "\n";

            if (i < ibl) i += 7;
            else if (i < ine) i += 5;
            else if (i < ifs) i += 1;
            else i += 5;
        }
        f << setw(width) << s.faceLoc(s.nVolumes()) * 100.0 << "\n";
        
        f.close();
        
    }*/
};

/**
 * Generic post processor.
 */
template <typename Process>
class PostProcessor
{
public:
    static void process(StagLineSol& s) {
        Process p;
        p.preLoop(s.getMixAtCellPoint(0));
        
        for (int i_cell = 0; i_cell < s.nVolumes(); ++i_cell)
            p.loop(s.cellLoc(i_cell)-s.faceLoc(0), s.getMixAtCellPoint(i_cell));
        
        p.postLoop();
    }
};

template <typename T>
class SpeciesProcessor
{
public:
    void preLoop(Mutation::Mixture& mix) {
        f.open(static_cast<T&>(*this).name().c_str());
        f << setw(15) << "x";
        for (int j = 0; j < mix.nSpecies(); ++j)
            f << setw(15) << mix.speciesName(j);
        f << "\n";
    }
    
    void loop(double x, Mutation::Mixture& mix) {
        f << setw(15) << x;
        for (int i = 0; i < mix.nSpecies(); ++i)
            f << setw(15) << static_cast<T&>(*this).value(mix, i);
        f << "\n";
    }
    
    void postLoop() {
        f.close();
    }
    
private:
    std::ofstream f;
};

/**
 * Stefan-Maxwell diffusion velocities.
 */
//class StefanMaxwell : public SpeciesProcessor<StefanMaxwell>
//{
//public:
//    std::string name() {
//        return "Vi_SM.dat";
//    }
//    double value(Mutation::Mixture& mix, int i) {
//        mix.stefanMaxwell
//    }
//};

/**
 * Creates rop.dat with reaction rates of progress
 */
class RatesOfProgress
{
public:
    void preLoop(Mutation::Mixture& mix) {
        f.open("rop.dat");
        f << setw(15) << "x";
        for (int j = 0; j < mix.nReactions(); ++j)
            f << setw(20) << '"' + mix.reactions()[j].formula() + '"';
        f << "\n";
    }

    void loop(double x, Mutation::Mixture& mix) {
        static std::vector<double> rop(mix.nReactions());
        mix.netRatesOfProgress(&rop[0]);
        f << setw(15) << x;
        for (int i = 0; i < mix.nReactions(); ++i)
            f << setw(20) << rop[i];
        f << "\n";
    }

    void postLoop() {
        f.close();
    }

private:
    std::ofstream f;
};

/**
 * Creates omegai.dat with species production rates due to reactions.
 */
class ProductionRates : public SpeciesProcessor<ProductionRates>
{
public:
    std::string name() {
        return "omegai.dat";
    }
    double value(Mutation::Mixture& mix, int i) {
        static std::vector<double> wdot(mix.nSpecies());
        mix.netProductionRates(&wdot[0]);
        return wdot[i];
    }
};



/**
 * Creates ni.dat with species number densities.
 */
class NumberDensities : public SpeciesProcessor<NumberDensities>
{
public:
    std::string name() {
        return "ni.dat";
    }
    double value(Mutation::Mixture& mix, int i) { 
        return mix.numberDensity()*std::max(mix.X()[i], 1.0e-99);
    }
};

/**
 * Creates rhoi.dat with species number densities.
 */
class Densities : public SpeciesProcessor<Densities>
{
public:
    std::string name() {
        return "rhoi.dat";
    }
    double value(Mutation::Mixture& mix, int i) { 
        return mix.density()*std::max(mix.Y()[i], 1.0e-99);
    }
};

/**
 * Creates xi.dat with species mole fractions.
 */
class MoleFractions : public SpeciesProcessor<MoleFractions>
{
public:
    std::string name() {
        return "xi.dat";
    }
    double value(Mutation::Mixture& mix, int i) { 
        return std::max(mix.X()[i], 1.0e-99);
    }
};

/**
 * Creates yi.dat with species mass fractions.
 */
class MassFractions : public SpeciesProcessor<MassFractions>
{
public:
    std::string name() {
        return "yi.dat";
    }
    double value(Mutation::Mixture& mix, int i) { 
        return std::max(mix.Y()[i], 1.0e-99);
    }
};

class TransferSourceTerm
{
public:

    void loop(double x, Mutation::Mixture& mix) {
        f << setw(15) << x;
        f << setw(20) << mp_omega->source();
        f << "\n";
    }

    void postLoop() {
        f.close();
    }

protected:

    void _preLoop(Mutation::Mixture& mix, const std::string& name) {
        f.open((name+".dat").c_str());
        f << setw(15) << "x";

        mp_omega =
            Mutation::Utilities::Config::Factory<
                Mutation::Transfer::TransferModel>::create(
                name, mix);
    }

private:

    std::ofstream f;
    Mutation::Transfer::TransferModel* mp_omega;
};

class OmegaVT : public TransferSourceTerm
{
public:
    void preLoop(Mutation::Mixture& mix) {
        _preLoop(mix, "OmegaVT");
    }
};

class OmegaCV : public TransferSourceTerm
{
public:
    void preLoop(Mutation::Mixture& mix) {
        _preLoop(mix, "OmegaCV");
    }
};

class OmegaET : public TransferSourceTerm
{
public:
    void preLoop(Mutation::Mixture& mix) {
        _preLoop(mix, "OmegaET");
    }
};

class OmegaCE : public TransferSourceTerm
{
public:
    void preLoop(Mutation::Mixture& mix) {
        _preLoop(mix, "OmegaCE");
    }
};

class OmegaCElec : public TransferSourceTerm
{
public:
    void preLoop(Mutation::Mixture& mix) {
        _preLoop(mix, "OmegaCElec");
    }
};

class OmegaI : public TransferSourceTerm
{
public:
    void preLoop(Mutation::Mixture& mix) {
        _preLoop(mix, "OmegaI");
    }
};



  
/**
 * Entry point for the program.  Loads solution and runs post processors.
 */
int main(int argc, char* argv[])
{
    // Load the solution
    StagLineSol solution(argc > 1 ? argv[1] : "input");
    
    // Run post processors on the solution
    PostProcessor<NumberDensities>::process(solution);
    PostProcessor<Densities>::process(solution);
    PostProcessor<MoleFractions>::process(solution);
    PostProcessor<MassFractions>::process(solution);
    PostProcessor<ProductionRates>::process(solution);
    PostProcessor<RatesOfProgress>::process(solution);
    PostProcessor<OmegaVT>::process(solution);
    PostProcessor<OmegaCV>::process(solution);
    PostProcessor<OmegaET>::process(solution);
    PostProcessor<OmegaCE>::process(solution);
    PostProcessor<OmegaCElec>::process(solution);
    PostProcessor<OmegaI>::process(solution);
    

    RadiationProcessor::process(solution);
    return 0;
}
