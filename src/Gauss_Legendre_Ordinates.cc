#include "Gauss_Legendre_Ordinates.hh"

#include "Gauss_Legendre.hh"
#include "Legendre_Polynomial.hh"

Gauss_Legendre_Ordinates::
Gauss_Legendre_Ordinates(int dimension,
                         int number_of_moments,
                         int number_of_ordinates): 
    dimension_(dimension),
    number_of_moments_(number_of_moments),
    number_of_ordinates_(number_of_ordinates)
{
    gauss_legendre_vec(number_of_ordinates, ordinates_, weights_);
}

double Gauss_Legendre_Ordinates::
moment(int mom,
       int ord)
{
    return legendre_polynomial(mom, ordinates_[ord]);
}
