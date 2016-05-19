#ifndef Finite_Element_hh
#define Finite_Element_hh

#include <vector>
#include "pugixml.hpp"

using std::vector;

/*
  A single finite element
*/
class Finite_Element
{
public:

    // Constructor
    Finite_Element(int dimension,
                   int number_of_nodes,
                   vector<double> const &node_positions);

    // Number of dimensions
    int dimension()
    {
        return dimension_;
    }

    // Number of nodes
    int number_of_nodes()
    {
        return number_of_nodes_;
    }

    // Length of element in each dimension
    vector<double> const &length() const
    {
        return length_;
    }

    // Positions of individual nodes
    vector<double> const &node_positions() const
    {
        return node_positions_;
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
    void output(pugi::xml_node &element_node) const;

private:

    int dimension_;
    int number_of_nodes_;
    vector<double> length_;
    vector<double> node_positions_;
};

#endif
