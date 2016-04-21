#include "Finite_Element.hh"

#include <vector>

#include "Check.hh"
#include "XML_Child_Value.hh"

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

void Finite_Element::
output(pugi::xml_node &element_node) const
{
    append_child(element_node, dimension_, "dimension");
    append_child(element_node, number_of_nodes_, "number_of_nodes"); 
    append_child(element_node, length_, "length", "dimension");
    append_child(element_node, node_positions_, "node_positions", "node");
}
