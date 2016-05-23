#ifndef Spatial_Discretization_hh
#define Spatial_Discretization_hh

#include <vector>
#include "pugixml.hpp"

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
    enum Geometry
    {
        SLAB,
        SPHERE,
        CYLINDER,
        RECTANGLE,
        CUBOID
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

    // Material number for each cell
    virtual vector<int> const &material() const = 0;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const = 0;
    
protected:
    
    int dimension_;
    Geometry geometry_;
};

#endif
