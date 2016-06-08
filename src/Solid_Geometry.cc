#include "Solid_Geometry.hh"

#include <numeric_limits>

using namespace std;

Solid_Geometry::
Solid_Geometry(int dimension,
               vector<shared_ptr<Surface> > const &surfaces,
               vector<shared_ptr<Region> > const &regions):
    surfaces_(surfaces),
    regions_(regions)
{
}

int Solid_Geometry::
find_region(vector<double> const &particle_position)
{
    Region::Relation relation;
    
    for (int i = 0; i < regions_.size(); ++i)
    {
        relation = regions_[i]->relation(particle_position);
        
        if (relation == Region::Relation::INSIDE)
        {
            return i;
        }
    }
    
    return -1; // no region found: outside of problem
}

int Solid_Geometry::
next_intersection(vector<double> const &particle_position,
                  vector<double> const &particle_direction,
                  double &distance,
                  vector<double> &position) const
{
    int best_surface = -1;
    distance = numeric_limits<double>::max();
    position.resize(dimension_);

    for (int i = 0; i < surfaces_.size(); ++i)
    {
        double current_distance;
        vector<double> current_position;
        
        if(intersection(particle_position,
                        particle_direction,
                        current_distance,
                        current_position))
        {
            if (current_distance < best_distance)
            {
                best_surface = i;
                distance = current_distance;
                position = current_position;
            }
        }
    }
    
    return best_surface;
}
