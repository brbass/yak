#ifndef Spatial_Discretization_hh
#define Spatial_Discretization_hh

class Spatial_Discretization
{
public:
    
    enum Geometry
    {
        SLAB,
        SPHERICAL,
        CYLINDRICAL
    };
    
    virtual int number_of_points() = 0;
    virtual int number_of_cells() = 0;
    virtual int number_of_nodes() = 0;
    virtual int dimension() = 0;
};

#endif
