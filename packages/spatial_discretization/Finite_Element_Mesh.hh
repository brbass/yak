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
    enum class Element_Type
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
    virtual int number_of_points() override
    {
        return number_of_points_;
    }

    // Number of spatial elements
    virtual int number_of_cells() override
    {
        return number_of_elements_;
    }

    virtual int number_of_materials() override
    {
        return number_of_materials_;
    }
    
    // Number of elements
    virtual int number_of_elements()
    {
        return number_of_elements_;
    }

    // Number of elements on the problem boundary
    virtual int number_of_boundary_cells() override
    {
        return number_of_boundary_elements_;
    }

    // Global number of nodes on the boundary
    virtual int number_of_boundary_points() override
    {
        return number_of_boundary_elements_;
    }

    virtual int number_of_internal_points() override
    {
        return number_of_internal_elements_;
    }
    
    // Number of nodes per element
    virtual int number_of_nodes() override
    {
        return number_of_nodes_;
    }

    // Which nodes are on the boundary for the boundary cells
    virtual vector<bool> const &boundary_nodes() const override
    {
        return boundary_nodes_;
    }

    // List of elements on the boundary
    virtual vector<int> const &boundary_cells() const override
    {
        return boundary_elements_;
    }

    // List of elements not on the boundary
    virtual vector<int> const &internal_cells() const override
    {
        return internal_elements_;
    }

    // Material number of each element
    virtual vector<int> const &material() const override
    {
        return material_;
    }

    virtual vector<double> const &boundary_normal() const override
    {
        return surface_normal_;
    }
    
    // Type of finite element
    virtual Element_Type element_type()
    {
        return element_type_;
    }

    // Return specific element
    virtual Finite_Element const &elements(int element) const
    {
        return elements_[element];
    }

    virtual vector<Cell_Type> const &cell_type() const override
    {
        return cell_type_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;
    
private:

    int number_of_elements_;
    int number_of_materials_;
    int number_of_nodes_;
    int number_of_points_;
    int number_of_internal_elements_;
    int number_of_boundary_elements_;
    
    vector<bool> boundary_nodes_;
    vector<int> boundary_elements_;
    vector<int> internal_elements_;
    vector<int> material_;
    vector<double> node_positions_;
    vector<double> surface_normal_;
    
    Element_Type element_type_;
    
    vector<Finite_Element> elements_;
    vector<Cell_Type> cell_type_;
};

#endif
