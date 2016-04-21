#include "Gauss_Legendre_Ordinates.hh"

#include "Check.hh"
#include "Gauss_Legendre.hh"
#include "XML_Child_Value.hh"

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

void Gauss_Legendre_Ordinates::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node gauss = output_node.append_child("gauss_legendre_ordinates");
    
    append_child(gauss, dimension_, "dimension");
    append_child(gauss, number_of_moments_, "number_of_moments");
    append_child(gauss, number_of_ordinates_, "number_of_ordinates");
    append_child(gauss, ordinates_, "ordinates", "ordinate");
    append_child(gauss, weights_, "weights", "ordinate");
}

