#include "Finite_Element.hh"

#include <vector>

using namespace std;

Finite_Element::
Finite_Element(int number_of_nodes,
               vector<double> &node_positions):
    number_of_nodes_(number_of_nodes),
    node_positions_(node_positions)
{
}
