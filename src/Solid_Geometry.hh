#ifndef Solid_Geometry_hh
#define Solid_Geometry_hh

#include "Region.hh"
#include "Surface.hh"

/* 
   Solid geometry for Monte Carlo and RBF calculations
*/
class Solid_Geometry
{
public:

    Solid_Geometry(int dimension,
                   vector<shared_ptr<Surface> > const &surfaces,
                   vector<shared_ptr<Region> > const &regions);

    shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    shared_ptr<Region> const &region(int r) const
    {
        return regions_[r];
    }

    int find_region(vector<double> const &particle_position);
    
    int next_intersection(vector<double> const &particle_position,
                          vector<double> const &particle_direction,
                          double &distance,
                          vector<double> &position) const;
    
private:

    int dimension_;
    vector<shared_ptr<Surface> > surfaces_;
    vector<shared_ptr<Region> > regions_;
};

#endif
