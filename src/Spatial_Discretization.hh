#ifndef Spatial_Discretization_hh
#define Spatial_Discretization_hh

#include <vector>

using std::vector;

class Spatial_Discretization
{
public:

    enum Geometry
    {
        SLAB,
        SPHERICAL,
        CYLINDRICAL
    };
    
    Spatial_Discretization(int dimension,
                           Geometry geometry);
    
    virtual int number_of_points() = 0;
    virtual int number_of_boundary_points() = 0;
    virtual int number_of_cells() = 0;
    virtual int number_of_boundary_cells() = 0;
    virtual int number_of_nodes() = 0;
    virtual int dimension()
    {
        return dimension_;
    }
    virtual Geometry geometry()
    {
        return geometry_;
    }
    virtual vector<int> const &boundary_cells() const = 0;
    virtual vector<int> const &internal_cells() const = 0;
    virtual vector<int> const &material() const = 0;
    
protected:
    
    int dimension_;
    Geometry geometry_;
};

#endif
