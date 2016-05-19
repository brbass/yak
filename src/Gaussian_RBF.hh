#ifndef Gaussian_RBF_hh
#define Gaussian_RBF_hh

#include "RBF.hh"

/*
  Radial basis function of the form exp(-c^2 r^2)
*/
class Gaussian_RBF : public RBF
{
public:

    // Constructor
    Gaussian_RBF(int number_of_dimensions,
                 vector<double> const &position,
                 vector<double> const &shape_parameter);

    // Value of basis function at the point r
    virtual double basis(vector<double> const &r) const;
    
    // Derivative of basis function at the point r
    virtual double dbasis(int dim,
                          vector<double> const &r) const;

    // Second derivative of the basis function at the point r
    virtual double ddbasis(int dim,
                           vector<double> const &r) const;
    
};

#endif
