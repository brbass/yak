#include "Lebedev_Quadrature.hh"

#include "Check.hh"
#include "Quadrature_Rule.hh"
#include "XML_Child_Value.hh"

Lebedev_Quadrature::
Lebedev_Quadrature(int dimension,
                   int number_of_moments,
                   int rule):
    Angular_Discretization(dimension,
                           number_of_moments,
                           Quadrature_Rule::lebedev(rule,
                                                    dimension,
                                                    ordinates_,
                                                    weights_))
{
    Quadrature_Rule::lebedev(rule,
                             dimension,
                             ordinates_,
                             weights_);
    
    check_class_invariants();
}

void Lebedev_Quadrature::
check_class_invariants() const
{
    Assert(ordinates_.size() == number_of_ordinates_ * dimension_);
    Assert(weights_.size() == number_of_ordinates_);
}

void Lebedev_Quadrature::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node gauss = output_node.append_child("gauss_legendre_ordinates");
    
    append_child(gauss, dimension_, "dimension");
    append_child(gauss, number_of_moments_, "number_of_moments");
    append_child(gauss, number_of_ordinates_, "number_of_ordinates");
    append_child(gauss, ordinates_, "ordinates", "dimension-ordinate");
    append_child(gauss, weights_, "weights", "ordinate");
}

