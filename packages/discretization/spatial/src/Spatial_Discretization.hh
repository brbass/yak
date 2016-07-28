#ifndef Spatial_Discretization_hh
#define Spatial_Discretization_hh

#include <vector>
#include "pugixml.hh"

using std::vector;

/*
  Pure virtual class for a spatial discretization

  Cells are the unit of physical data, such as cross sections
  Points are the solution points for the data
  Nodes are the points within each cell
*/
class Spatial_Discretization
{
public:

    // Geometry of problem
    enum class Geometry
    {
        SLAB,
        SPHERE,
        CYLINDER,
        CARTESIAN
    };

    // Constructor
    Spatial_Discretization(int dimension,
                           Geometry geometry);

    // Number of spatial points to solve for in the problem
    virtual int number_of_points() = 0;

    // Number of spatial points not on the boundary
    virtual int number_of_internal_points() = 0;
    
    // Number of spatial points on the boundary
    virtual int number_of_boundary_points() = 0;

    // Number of cells transitioning from one material to another
    virtual int number_of_transition_points()
    {
        return 0;
    }
    
    // Number of cells
    virtual int number_of_cells() = 0;

    // Number of cells on the boundary
    virtual int number_of_boundary_cells() = 0;

    // Number of nodes per cell
    virtual int number_of_nodes() = 0;

    // Number of spatial dimensions
    virtual int dimension()
    {
        return dimension_;
    }
    
    virtual int number_of_materials() = 0;
    
    // Geometry of problem
    virtual Geometry geometry()
    {
        return geometry_;
    }

    // Which nodes in each boundary cell are themselves on the boundary
    virtual vector<bool> const &boundary_nodes() const = 0;

    // Spatial points on the problem boundary
    virtual vector<int> const &boundary_cells() const = 0;

    // Spatial points not on the problem boundary
    virtual vector<int> const &internal_cells() const = 0;

    // Spatial points transitioning from one material to another
    virtual vector<int> const &transition_cells() const
    {
        return empty_int_;
    }
    
    // Material number for each cell
    virtual vector<int> const &material() const = 0;

    // Surface normal values
    virtual vector<double> const &boundary_normal() const = 0;

    virtual vector<double> const &transition_normal() const
    {
        return empty_double_;
    }
    
    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const = 0;
    
protected:
    
    int dimension_;
    Geometry geometry_;

    vector<int> empty_int_;
    vector<double> empty_double_;
};

#endif
