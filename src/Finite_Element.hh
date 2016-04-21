#ifndef Finite_Element_hh
#define Finite_Element_hh

#include <vector>
#include "pugixml.hpp"

using std::vector;

class Finite_Element
{
public:

    Finite_Element(int dimension,
                   int number_of_nodes,
                   vector<double> const &node_positions);

    int dimension()
    {
        return dimension_;
    }
    int number_of_nodes()
    {
        return number_of_nodes_;
    }
    vector<double> const &length() const
    {
        return length_;
    }
    vector<double> const &node_positions() const
    {
        return node_positions_;
    }

    void check_class_invariants() const;
    void output(pugi::xml_node &element_node) const;

private:

    int dimension_;
    int number_of_nodes_;
    vector<double> length_;
    vector<double> node_positions_;
};

#endif
