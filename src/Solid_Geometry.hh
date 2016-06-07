#ifndef Solid_Geometry_hh
#define Solid_Geometry_hh

#include "Region.hh"
#include "Surface.hh"

class Solid_Geometry
{
public:

    Solid_Geometry();

private:

    vector<Surface> surfaces_;
    vector<Region> regions_;
    vector<int> material_;
    
};

#endif
