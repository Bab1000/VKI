
#include "mutation++.h"
#include <string>

/**
 * Convenience class which represents a single flowfield point.
 */
class StagLinePoint
{
public:
    StagLinePoint(const Mutation::Mixture& mix);
    
    // Getter functions
    double* const rhoi() { return &m_data[0]; }
    double* const Tk() { return &m_data[m_ns+2]; }
    double& u() { return m_data[m_ns]; }
    double& v() { return m_data[m_ns+1]; }
    
    /**
     * Sets the mixture state to this point.
     */
    Mutation::Mixture& setMixtureState(Mutation::Mixture& mix) {
        mix.setState(rhoi(), Tk(), 1);
        return mix;
    }
    
    // Some helper math functions
    StagLinePoint operator+ (const StagLinePoint& p) {
        StagLinePoint point(*this);
        for (int i = 0; i < m_ns+m_nt+2; ++i)
            point.m_data[i] += p.m_data[i];
        return point;
    }
    
    StagLinePoint operator* (double a) {
        StagLinePoint point(*this);
        for (int i = 0; i < m_ns+m_nt+2; ++i)
            point.m_data[i] *= a;
        return point;
    }
    
private:
    int m_ns, m_nt;
    std::vector<double> m_data;
};

class StagLineSol
{
public:
    StagLineSol(const std::string& filename);
    ~StagLineSol();
    
    size_t nVolumes() const { return m_centers.size(); }
    size_t nFaces()   const { return m_faces.size(); }
    size_t nSpecies() const { return mp_mix->nSpecies(); }
    size_t nEnergyEqns() const { return mp_mix->nEnergyEqns(); }
    
    double cellLoc(int i) const { return m_centers[i]; }
    double faceLoc(int i) const { return m_faces[i]; }
    
    Mutation::Mixture& getMixAtFacePoint(int i) {
        return (m_face_data[i].setMixtureState(*mp_mix));
    }
    Mutation::Mixture& getMixAtCellPoint(int i) {
        return (m_cell_data[i].setMixtureState(*mp_mix));
    }
    
private:
    /**
     * Reads the input file to figure the solution and mixture names.
     * @param filename - name of the input file
     */
    void parseInputFile(const std::string& filename);
    
    /**
     * Gets the mesh information.
     */
    void readMesh();
    
    /**
     * Gets the solution from the restart file.
     */
    void readRestartFile();

private:
    
    std::string m_case_name;
    std::string m_mixture_name;
    std::string m_state_model;
    
    Mutation::Mixture* mp_mix;
    
    std::vector<double> m_faces;
    std::vector<double> m_centers;
    
    std::vector<StagLinePoint> m_cell_data;
    std::vector<StagLinePoint> m_face_data;
};
