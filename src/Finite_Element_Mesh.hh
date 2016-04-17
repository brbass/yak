#ifndef Finite_Element_Mesh_hh
#define Finite_Element_Mesh_hh

#include <string>
#include <vector>

#include "Finite_Element.hh"
#include "Spatial_Discretization.hh"

using std::string;
using std::vector;

class Finite_Element_Mesh : public Spatial_Discretization
{
public:
    
    enum Finite_Element_Type
    {
        DFEM,
        CFEM
    };

    Finite_Element_Mesh(int dimension,
                        int number_of_elements,
                        int number_of_nodes,
                        Geometry geometry,
                        Finite_Element_Type element_type,
                        vector<double> const &node_positions);
    
    virtual int number_of_points()
    {
        return number_of_points_;
    }
    virtual int number_of_cells()
    {
        return number_of_elements_;
    }
    virtual int number_of_nodes()
    {
        return number_of_nodes_;
    }
    Finite_Element_Type finite_element_type()
    {
        return element_type_;
    }
    Finite_Element const &elements(int element) const
    {
        return elements_[element];
    }
    
    void check_class_invariants() const;
    
private:

    int number_of_elements_;
    int number_of_nodes_;
    int number_of_points_;
    
    Finite_Element_Type element_type_;
    
    vector<Finite_Element> elements_;
};

#endif
