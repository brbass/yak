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
    
    virtual int dimension()
    {
        return dimension_;
    }
    virtual int number_of_moments()
    {
        return number_of_moments_;
    }
    virtual int number_of_ordinates()
    {
        return number_of_ordinates_;
    }
    virtual vector<double> const &ordinates()
    {
        return ordinates_;
    }
    virtual vector<double> const &weights()
    {
        return weights_;
    }
    virtual double moment(int mom,
                          int ord);

    void check_class_invariants();
    
private:

    int dimension_;
    int number_of_ordinates_;
    int number_of_moments_;

    vector<double> ordinates_;
    vector<double> weights_;
};

#endif
