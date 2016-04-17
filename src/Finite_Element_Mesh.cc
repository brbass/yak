#include "Finite_Element_Mesh.hh"

#include <string>
#include <vector>

#include "Check.hh"

using namespace std;

Finite_Element_Mesh::
Finite_Element_Mesh(int dimension, 
                    int number_of_elements,
                    int number_of_nodes,
                    Geometry geometry,
                    Finite_Element_Type element_type,
                    vector<double> const &node_positions):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_elements_(number_of_elements),
    number_of_nodes_(number_of_nodes),
    number_of_boundary_elements_(2),
    element_type_(element_type)
{

    switch(element_type_)
    {
    case CFEM:
        number_of_points_ = number_of_elements_ + 1;
        break;
    case DFEM:
        number_of_points_ = number_of_elements_ * number_of_nodes_;
    }

    for (int i = 0; i < number_of_elements_; ++i)
    {
        int k1 = dimension * number_of_nodes_ * i;

        vector<double> local_node_positions(number_of_nodes * dimension);
        for (int d = 0; d < dimension_; ++d)
        {
            for (int n = 0; n < number_of_nodes_; ++n)
            {
                int k1 = d + dimension * (n + number_of_nodes_ * i);
                int k2 = d + dimension * n;

                local_node_positions[k2] = node_positions[k1];
            }
        }
        
        elements_.emplace_back(dimension,
                               number_of_nodes,
                               local_node_positions);
    }

    boundary_elements_.push_back(0);
    boundary_elements_.push_back(number_of_elements_ - 1);
    
    for (int i = 1; i < number_of_elements_ - 1; ++i)
    {
        internal_elements_.push_back(i);
    }
    
    check_class_invariants();
}

void Finite_Element_Mesh::
check_class_invariants() const
{
    Assert(elements_.size() == number_of_elements_);
    Assert(boundary_elements_.size() == number_of_boundary_elements_);
    Assert(internal_elements_.size() == number_of_elements_ - number_of_boundary_elements_);
}
