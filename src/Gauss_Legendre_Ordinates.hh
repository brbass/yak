#ifndef Gauss_Legendre_Ordinates_hh
#define Gauss_Legendre_Ordinates_hh

#include <vector>

#include "Angular_Discretization.hh"

using std::vector;

class Gauss_Legendre_Ordinates : public Angular_Discretization
{
public:
    
    Gauss_Legendre_Ordinates(int dimension,
                             int number_of_moments,
                             int number_of_ordinates);
    
    virtual vector<double> const &ordinates() const
    {
        return ordinates_;
    }
    virtual vector<double> const &weights() const
    {
        return weights_;
    }
    
    void check_class_invariants() const;
    
private:

    vector<double> ordinates_;
    vector<double> weights_;
};

#endif
