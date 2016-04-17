#include "Angular_Discretization.hh"

#include <cmath>
#include <vector>

#include "Check.hh"
#include "Legendre_Polynomial.hh"

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

double Angular_Discretization::
moment(int mom,
       int ord)
{
    switch(dimension_)
    {
    case 1:
    {
        double mu = ordinates()[ord];
        return legendre_polynomial(mom, mu);
    }
    default:
        AssertMsg(false, "Only one dimension supported");
        return 0;
    }
}
