#include "Angular_Discretization.hh"
#include "Legendre_Polynomial.hh"

#include <vector>

using namespace std;

Angular_Discretization::
Angular_Discretization(int dimension, 
                       int number_of_moments,
                       int number_of_ordinates):
    dimension_(dimension),
    number_of_moments_(number_of_moments),
    number_of_ordinates_(number_of_ordinates)
{
    switch(dimension_)
    {
    case 1:
        angular_normalization_ = 2;
        break;
    default:
        angular_normalization_ = 4 * M_PI;
        break;
    }
}

virtual double Angular_Discretization::
moment(int mom,
       int ord)
{
    switch(dimension_)
    {
    case 1:
        return legendre_polynomial(mom, ordinates()[ord]);
    default:
    }
}
