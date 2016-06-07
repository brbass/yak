#ifndef Solid_Geometry_hh
#define Solid_Geometry_hh

#include "Region.hh"
#include "Surface.hh"

class Solid_Geometry
{
public:

    Solid_Geometry();

    shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    shared_ptr<Region> const &region(int r) const
    {
        return regions_[r];
    }
    
private:

    vector<shared_ptr<Surface> > surfaces_;
    vector<shared_ptr<Region> > regions_;
};

#endif
