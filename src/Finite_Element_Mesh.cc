#include "Finite_Element_Mesh.hh"

#include <string>
#include <vector>

#include "Check.hh"

using namespace std;

Finite_Element_Mesh::
Finite_Element_Mesh(int dimension, 
                    int number_of_elements,
                    int number_of_nodes,
                    string element_type,
                    vector<double> &node_positions):
    dimension_(dimension),
    number_of_elements_(number_of_elements),
    number_of_nodes_(number_of_nodes)
{
    set_element_type(element_type);

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
        
        elements_.emplace_back(number_of_nodes,
                               local_node_positions);
    }

    check_class_invariants();
}

Finite_Element_Mesh::
Finite_Element_Mesh(int dimension,
                    int number_of_elements,
                    int number_of_nodes,
                    double length,
                    string element_type):
    dimension_(dimension),
    number_of_elements_(number_of_elements),
    number_of_nodes_(number_of_nodes)    
{
    Assert(dimension == 1);

    set_element_type(element_type);

    double dx = length / number_of_elements;
    double ddx = dx / (number_of_nodes - 1);
    
    for (int i = 0; i < number_of_elements_; ++i)
    {
        vector<double> node_positions(number_of_nodes);
        
        for (int n = 0; n < number_of_nodes_; ++n)
        {
            double position = dx * i + n * ddx;
            
            node_positions[n] = position;
        }
        
        elements_.emplace_back(number_of_nodes,
                               node_positions);
    }

    check_class_invariants();
}

void Finite_Element_Mesh::
set_element_type(string element_type)
{
    if (element_type == "continuous")
    {
        element_type_ = CFEM;
    }
    else 
    {
        element_type_ = DFEM;
    }
}

void Finite_Element_Mesh::
check_class_invariants()
{
    Assert(elements_.size() == number_of_elements_);
}
