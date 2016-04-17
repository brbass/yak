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

    Finite_Element_Mesh(int dimension,
                        int number_of_elements,
                        int number_of_nodes,
                        string element_type,
                        vector<double> const &node_positions);
    
    Finite_Element_Mesh(int dimension,
                        int number_of_elements,
                        int number_of_nodes,
                        double length, 
                        string element_type);
    
    enum Finite_Element_Type
    {
        DFEM,
        CFEM
    };
    
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
    virtual int dimension()
    {
        return dimension_;
    }
    Finite_Element const &elements(int element) const
    {
        return elements_[element];
    }

    void check_class_invariants() const;
    
private:
    
    void set_element_type(string element_type);

    int dimension_;
    int number_of_elements_;
    int number_of_nodes_;
    int number_of_points_;
    
    Finite_Element_Type element_type_;

    vector<Finite_Element> elements_;
};

#endif
