#include "Gaussian_RBF.hh"

#include <cmath>

Gaussian_RBF::
Gaussian_RBF(int number_of_dimensions,
             vector<double> const &position,
             vector<double> const &shape_parameter):
    RBF(number_of_dimensions,
        position,
        shape_parameter)
{
}

double Gaussian_RBF::
basis(vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    
    return exp(-dist2);
}

double Gaussian_RBF::
dbasis(int dim,
       vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double kr2 = shape_[dim] * shape_[dim] * (r[dim] - position_[dim]);
    
    return -2 * kr2 * exp(-dist2);
}

double Gaussian_RBF::
ddbasis(int dim,
        vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double distdim = r[dim] - position_[dim];
    double k2r2 = shape_[dim] * shape_[dim] * distdim * distdim;
    
    return 2 * (2 * k2r2 - 1) * shape_[dim] * shape_[dim] * exp(-dist2);
}
