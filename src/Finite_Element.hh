#ifndef Finite_Element_hh
#define Finite_Element_hh

#include <vector>

using std::vector;

class Finite_Element
{
public:

    Finite_Element(int number_of_nodes,
                   vector<double> const &node_positions);

    int number_of_nodes()
    {
        return number_of_nodes_;
    }
    vector<double> const &node_positions() const
    {
        return node_positions_;
    }

    void check_class_invariants() const;

private:

    int number_of_nodes_;
    vector<double> node_positions_;
};

#endif
