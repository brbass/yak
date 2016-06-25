#include "Gauss_Legendre_Quadrature.hh"

#include "Check.hh"
#include "Quadrature_Rule.hh"
#include "XML_Functions.hh"

Gauss_Legendre_Quadrature::
Gauss_Legendre_Quadrature(int dimension,
                          int number_of_moments,
                          int number_of_ordinates):
    Angular_Discretization(dimension,
                           number_of_moments,
                           number_of_ordinates)
{
    Quadrature_Rule::gauss_legendre(number_of_ordinates, ordinates_, weights_);
    
    check_class_invariants();
}

void Gauss_Legendre_Quadrature::
check_class_invariants() const
{
    Assert(ordinates_.size() == number_of_ordinates_);
    Assert(weights_.size() == number_of_ordinates_);
}

void Gauss_Legendre_Quadrature::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node gauss = output_node.append_child("gauss_legendre_ordinates");
    
    XML_Functions::append_child(gauss, dimension_, "dimension");
    XML_Functions::append_child(gauss, number_of_moments_, "number_of_moments");
    XML_Functions::append_child(gauss, number_of_ordinates_, "number_of_ordinates");
    XML_Functions::append_child(gauss, ordinates_, "ordinates", "ordinate");
    XML_Functions::append_child(gauss, weights_, "weights", "ordinate");
}

