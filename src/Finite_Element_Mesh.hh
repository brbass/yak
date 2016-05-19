#ifndef Finite_Element_Mesh_hh
#define Finite_Element_Mesh_hh

#include <string>
#include <vector>

#include "Finite_Element.hh"
#include "Spatial_Discretization.hh"

using std::string;
using std::vector;

/*
  Represents a general finite element mesh
  
  Only implemented in 1D for now
*/
class Finite_Element_Mesh : public Spatial_Discretization
{
public:

    // Discontinuous or continuous finite element
    enum Element_Type
    {
        DFEM,
        CFEM
    };

    // Constructor
    Finite_Element_Mesh(int dimension,
                        int number_of_elements,
                        int number_of_nodes,
                        Geometry geometry,
                        Element_Type element_type,
                        vector<int> const &material,
                        vector<double> const &node_positions);

    // Global number of nodes
    virtual int number_of_points()
    {
        return number_of_points_;
    }

    // Number of spatial elements
    virtual int number_of_cells()
    {
        return number_of_elements_;
    }

    // Number of elements
    virtual int number_of_elements()
    {
        return number_of_elements_;
    }

    // Number of elements on the problem boundary
    virtual int number_of_boundary_cells()
    {
        return number_of_boundary_elements_;
    }

    // Global number of nodes on the boundary
    virtual int number_of_boundary_points()
    {
        return number_of_boundary_elements_;
    }

    virtual int number_of_internal_points()
    {
        return number_of_internal_elements_;
    }
    
    // Number of nodes per element
    virtual int number_of_nodes()
    {
        return number_of_nodes_;
    }

    // Which nodes are on the boundary for the boundary cells
    virtual vector<bool> const &boundary_nodes() const
    {
        return boundary_nodes_;
    }

    // List of elements on the boundary
    virtual vector<int> const &boundary_cells() const
    {
        return boundary_elements_;
    }

    // List of elements not on the boundary
    virtual vector<int> const &internal_cells() const
    {
        return internal_elements_;
    }

    // Material number of each element
    virtual vector<int> const &material() const
    {
        return material_;
    }

    // Type of finite element
    Element_Type element_type()
    {
        return element_type_;
    }

    // Return specific element
    Finite_Element const &elements(int element) const
    {
        return elements_[element];
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const;
    
private:

    int number_of_elements_;
    int number_of_nodes_;
    int number_of_points_;
    int number_of_internal_elements_;
    int number_of_boundary_elements_;
    
    vector<bool> boundary_nodes_;
    vector<int> boundary_elements_;
    vector<int> internal_elements_;
    vector<int> material_;

    Element_Type element_type_;
    
    vector<Finite_Element> elements_;
};

#endif
