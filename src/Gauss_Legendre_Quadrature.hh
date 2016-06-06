#ifndef Gauss_Legendre_Quadrature_hh
#define Gauss_Legendre_Quadrature_hh

#include <vector>

#include "Angular_Discretization.hh"

using std::vector;

/*
  Holds Gauss Legendre quadrature
*/
class Gauss_Legendre_Quadrature : public Angular_Discretization
{
public:

    // Constructor
    Gauss_Legendre_Quadrature(int dimension,
                              int number_of_moments,
                              int number_of_ordinates);

    // Return Gauss-Legendre ordinates
    virtual vector<double> const &ordinates() const
    {
        return ordinates_;
    }
    
    // Return Gauss-Legendre weights
    virtual vector<double> const &weights() const
    {
        return weights_;
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const;

    virtual int reflect_ordinate(int o,
                                 vector<double> const &/*normal*/) const
    {
        return number_of_ordinates_ - o - 1;
    }
    
private:

    vector<double> ordinates_;
    vector<double> weights_;
};

#endif
