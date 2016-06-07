#include "Finite_Element_Mesh.hh"

#include <string>
#include <vector>

#include "Check.hh"
#include "XML_Functions.hh"

using namespace std;

Finite_Element_Mesh::
Finite_Element_Mesh(int dimension, 
                    int number_of_elements,
                    int number_of_nodes,
                    Geometry geometry,
                    Element_Type element_type,
                    vector<int> const &material,
                    vector<double> const &node_positions):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_elements_(number_of_elements),
    number_of_nodes_(number_of_nodes),
    number_of_internal_elements_(number_of_elements - 2),
    number_of_boundary_elements_(2),
    element_type_(element_type),
    material_(material),
    node_positions_(node_positions)
{
    Assert(dimension == 1);
    
    switch(element_type_)
    {
    case Element_Type::CFEM:
        number_of_points_ = number_of_elements_ + 1;
        break;
    case Element_Type::DFEM:
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
    
    boundary_nodes_.assign(number_of_nodes * number_of_boundary_elements_, false);
    boundary_nodes_[0 + number_of_nodes * 0] = true;
    boundary_nodes_[(number_of_nodes - 1) + number_of_nodes * (number_of_boundary_elements_ - 1)] = true;
    
    for (int i = 1; i < number_of_elements_ - 1; ++i)
    {
        internal_elements_.push_back(i);
    }

    surface_normal_.resize(2);
    surface_normal_[0] = -1;
    surface_normal_[1] = 1;
    
    check_class_invariants();
}

void Finite_Element_Mesh::
check_class_invariants() const
{
    Assert(elements_.size() == number_of_elements_);
    Assert(boundary_elements_.size() == number_of_boundary_elements_);
    Assert(boundary_nodes_.size() == number_of_boundary_elements_ * number_of_nodes_);
    Assert(internal_elements_.size() == number_of_elements_ - number_of_boundary_elements_);
    Assert(surface_normal_.size() == number_of_boundary_elements_ * dimension_);
}

void Finite_Element_Mesh::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node finite = output_node.append_child("finite_element_mesh");
    
    XML_Functions::append_child(finite, number_of_elements_, "number_of_elements");
    XML_Functions::append_child(finite, number_of_nodes_, "number_of_nodes");
    XML_Functions::append_child(finite, number_of_points_, "number_of_points");
    XML_Functions::append_child(finite, number_of_boundary_elements_, "number_of_boundary_elements");
    XML_Functions::append_child(finite, boundary_nodes_, "boundary_nodes", "element-node");
    XML_Functions::append_child(finite, boundary_elements_, "boundary_elements", "element");
    XML_Functions::append_child(finite, internal_elements_, "internal_elements", "element");
    XML_Functions::append_child(finite, material_, "material", "element");
    XML_Functions::append_child(finite, node_positions_, "node_positions", "element-node");
    
    for (int i = 0; i < number_of_elements_; ++i)
    {
        pugi::xml_node finite_element = finite.append_child("element");
        pugi::xml_attribute num = finite_element.append_attribute("number");
        num.set_value(to_string(i).c_str());
        
        elements_[i].output(finite_element);
    }
}
