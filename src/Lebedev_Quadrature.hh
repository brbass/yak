#ifndef Lebedev_Quadrature_hh
#define Lebedev_Quadrature_hh

#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"

using std::vector;

/*
  Holds Lebedev quadrature
*/
class Lebedev_Quadrature : public Angular_Discretization
{
public:

    // Constructor
    Lebedev_Quadrature(int dimension,
                       int number_of_moments,
                       int number_of_ordinates);

    // Return Lebedev ordinates
    virtual vector<double> const &ordinates() const
    {
        return ordinates_;
    }
    
    // Return Lebedev weights
    virtual vector<double> const &weights() const
    {
        return weights_;
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const;

    // Not yet implemented
    virtual int reflect_ordinate(int o,
                                 vector<double> const &n) const
    {
        AssertMsg(false, "not yet implemented");
        
        return o;
    }
    
private:

    vector<double> ordinates_;
    vector<double> weights_;
};

#endif
