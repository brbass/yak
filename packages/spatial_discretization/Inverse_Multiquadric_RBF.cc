#include "Inverse_Multiquadric_RBF.hh"

#include <cmath>

Inverse_Multiquadric_RBF::
Inverse_Multiquadric_RBF(int number_of_dimensions,
             vector<double> const &position,
             vector<double> const &shape_parameter):
    RBF(number_of_dimensions,
        position,
        shape_parameter)
{
}

double Inverse_Multiquadric_RBF::
basis(vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    
    return 1 / sqrt(1 + dist2);
}

double Inverse_Multiquadric_RBF::
dbasis(int dim,
       vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double kr2 = shape_[dim] * shape_[dim] * (r[dim] - position_[dim]);
    
    return -kr2 * pow(1 + dist2, -1.5);
}

double Inverse_Multiquadric_RBF::
ddbasis(int dim,
        vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double distdim = r[dim] - position_[dim];
    double k2r2 = shape_[dim] * shape_[dim] * distdim * distdim;
    
    return shape_[dim] * shape_[dim] * (1 + dist2 - 3 * k2r2) * pow(1 + dist2, -2.5);
}
