#ifndef LDFE_Quadrature_hh
#define LDFE_Quadrature_hh

#include "Angular_Discretization.hh"

/*
  Quadrature based on linear discontinuous finite elements
  
  Ref: Discrete-ordinates quadrature sets based on linear discontinuous finite elements
  Authors: Jarrell and Adams
*/
class LDFE_Quadrature : public Angular_Discretization
{
public:
    
    // Constructor
    LDFE_Quadrature(int dimension,
                    int number_of_scattering_moments,
                    int rule);
    
    // Return all ordinates
    virtual vector<double> const &ordinates() const
    {
        return ordinates_;
    }

    // Return weights
    virtual vector<double> const &weights() const
    {
        return weights_;
    }
    
    // X component of ordinates
    virtual vector<double> const &mu() const
    {
        return mu_;
    }

    // Y component of ordinates
    virtual vector<double> const &eta() const
    {
        return mu_;
    }
    
    // Z component of ordinates
    virtual vector<double> const &xi() const
    {
        return mu_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const;

private:

    void initialize_quadrature();
    void sort_quadrature();
    void initialize_1();
    void initialize_2();
    void initialize_3();
    
    int rule_;

    vector<double> mu_;
    vector<double> eta_;
    vector<double> xi_;
    vector<double> ordinates_;
    vector<double> weights_;
};

#endif
