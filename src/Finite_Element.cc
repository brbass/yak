#include "Finite_Element.hh"

#include <vector>

#include "Check.hh"

using namespace std;

Finite_Element::
Finite_Element(int dimension,
               int number_of_nodes,
               vector<double> const &node_positions):
    dimension_(dimension),
    number_of_nodes_(number_of_nodes),
    node_positions_(node_positions)
{
    Assert(dimension_ == 1);
    length_.push_back(node_positions_[number_of_nodes - 1] - node_positions[0]);
    
    check_class_invariants();
}

void Finite_Element::
check_class_invariants() const
{
    Assert(node_positions_.size() == number_of_nodes_);
}
