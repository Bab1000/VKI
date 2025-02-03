#include "StagLineSol.h"
#include <iostream>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using namespace Mutation;

//==============================================================================

StagLinePoint::StagLinePoint(const Mixture& mix)
{
    m_ns = mix.nSpecies();
    m_nt = mix.nEnergyEqns();
    m_data.assign(m_ns+m_nt+2, 0.0);
}

//==============================================================================

StagLineSol::StagLineSol(const std::string& filename)
{
    // First get the names of the mixture and solution
    parseInputFile(filename);
    cout << "Solution Name: " << m_case_name << endl;
    cout << "Mixture Name:  " << m_mixture_name << endl;
    
    // Load the mesh
    readMesh();
    cout << "# Volumes:     " << nVolumes() << endl;
    
    // Load the mixture
    MixtureOptions opts(m_mixture_name);
    opts.setStateModel(m_state_model);
    mp_mix = new Mixture(opts);
    
    // Load the solution from the restart file
    readRestartFile();
    
    // Interpolate the data to the cell faces (note the boundary faces are just
    // directly copied from the closest cell center for now)
    m_face_data.push_back(m_cell_data[0]);
    double v, v1, v2;
    for (int i = 0; i < nVolumes()-1; ++i) {
        v1 = 0.5*(m_faces[i]+m_faces[i+1]);
        v2 = 0.5*(m_faces[i+1]+m_faces[i+2]);
        v  = v1 + v2; v1 /= v; v2 /= v;
        
        m_face_data.push_back(m_cell_data[i]*v1+m_cell_data[i+1]*v2);
    }
    m_face_data.push_back(m_cell_data.back());
}

//==============================================================================

StagLineSol::~StagLineSol()
{
    if (mp_mix != NULL) delete mp_mix;
}

//==============================================================================

void StagLineSol::parseInputFile(const std::string& filename)
{
    std::ifstream f(filename.c_str());
    if (!f.is_open()) {
        cout << "Could not find input file: " << filename << endl;
        exit(1);
    }
    
    std::string line;
    getline(f, line);
    while (line != "") {
        if (line == "Simulation_Name") {
            getline(f, line);
            m_case_name = line;
        } else if (line == "Mixture") {
            getline(f, line);
            m_mixture_name = line;
        } else if (line == "State_Model") {
            getline(f, line);
            m_state_model = line;
        }
        getline(f, line);
    }
    
    f.close();
}

//==============================================================================

void StagLineSol::readMesh()
{
    // Open the mesh.dat file
    std::ifstream f("mesh.dat");
    if (!f.is_open()) {
        cout << "Could not find mesh.dat" << endl;
        exit(1);
    }
    
    // Read the faces
    int nfaces; f >> nfaces;
    m_faces.resize(nfaces);
    for (int i = 0; i < nfaces; ++i)
        f >> m_faces[i];
    
    f.close();
    
    // Get the cell centers
    m_centers.resize(nfaces-1);
    for (int i = 0; i < nfaces-1; ++i)
        m_centers[i] = 0.5*(m_faces[i] + m_faces[i+1]);
}

//==============================================================================

void StagLineSol::readRestartFile()
{
    // Open the restart file
    std::string filename = m_case_name+"_restartBC.dat";
    std::ifstream f(filename.c_str());
    if (!f.is_open()) {
        cout << "Could not find restart file: " << filename << endl;
        exit(1);
    }
    
    double P;
    for (int i = 0; i < nVolumes(); ++i) {
        m_cell_data.push_back(*mp_mix);
        
        // Species densities
        for (int k = 0; k < mp_mix->nSpecies(); ++k)
            f >> m_cell_data[i].rhoi()[k];
        
        // Velocities
        f >> m_cell_data[i].u();
        f >> m_cell_data[i].v();
        
        // Temperatures
        for (int k = 0; k < mp_mix->nEnergyEqns(); ++k)
            f >> m_cell_data[i].Tk()[k];
    }
    
    f.close();
}
