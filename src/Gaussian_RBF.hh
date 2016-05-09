#ifndef Gaussian_RBF_hh
#define Gaussian_RBF_hh

#include "RBF.hh"

class Gaussian_RBF : public RBF
{
public:

    Gaussian_RBF(int number_of_dimensions,
                 vector<double> const &position,
                 vector<double> const &shape_parameter);
    
    virtual double basis(vector<double> const &r) const;
    virtual double dbasis(int dim,
                          vector<double> const &r) const;
    virtual double ddbasis(int dim,
                           vector<double> const &r) const;
    
};

#endif
