#include "Gauss_Legendre_Ordinates.hh"

#include "Check.hh"
#include "Gauss_Legendre.hh"

Gauss_Legendre_Ordinates::
Gauss_Legendre_Ordinates(int dimension,
                         int number_of_moments,
                         int number_of_ordinates): 
    Angular_Discretization(dimension,
                           number_of_moments,
                           number_of_ordinates)
{
    gauss_legendre_vec(number_of_ordinates, ordinates_, weights_);
    
    check_class_invariants();
}

void Gauss_Legendre_Ordinates::
check_class_invariants() const
{
    Assert(ordinates_.size() == number_of_ordinates_);
    Assert(weights_.size() == number_of_ordinates_);
}
