#include "Finite_Element.hh"

#include <vector>

#include "Check.hh"

using namespace std;

Finite_Element::
Finite_Element(int number_of_nodes,
               vector<double> &node_positions):
    number_of_nodes_(number_of_nodes),
    node_positions_(node_positions)
{
    check_class_invariants();
}

void Finite_Element::
check_class_invariants()
{
    Assert(node_positions_.size() == number_of_nodes_);
}
