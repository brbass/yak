#include "Solid_Geometry.hh"

#include <limits>

using namespace std;

namespace // anonymous
{
    int get_number_of_materials(vector<shared_ptr<Region> > const &regions)
    {
        int number_of_regions = regions.size();

        int highest_material = 0;

        for (int r = 0; r < number_of_regions; ++r)
        {
            int region_material = regions[r]->material();

            if (region_material > highest_material)
            {
                highest_material = region_material;
            }
        }

        return highest_material + 1;
    }
} // end namespace anonymous

Solid_Geometry::
Solid_Geometry(int dimension,
               vector<shared_ptr<Surface> > const &surfaces,
               vector<shared_ptr<Region> > const &regions):
    dimension_(dimension),
    number_of_materials_(get_number_of_materials(regions)),
    surfaces_(surfaces),
    regions_(regions)
{
    delta_distance_ = 1000 * numeric_limits<double>::epsilon();
}

int Solid_Geometry::
find_region(vector<double> const &position) const
{
    Region::Relation relation;
    
    for (int i = 0; i < regions_.size(); ++i)
    {
        relation = regions_[i]->relation(position);
        
        if (relation == Region::Relation::INSIDE)
        {
            return i;
        }
    }
    
    return NO_REGION; // no region found: outside of problem
}

int Solid_Geometry::
find_surface(vector<double> const &position) const
{
    Surface::Relation relation;
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        relation = surfaces_[i]->relation(position);
        
        if (relation == Surface::Relation::EQUAL)
        {
            return i;
        }
    }
    
    return NO_SURFACE; // particle not on a surface
}

int Solid_Geometry::
next_geometric_intersection(vector<double> const &initial_position,
                            vector<double> const &initial_direction,
                            double &distance,
                            vector<double> &final_position) const
{
    int best_surface = NO_SURFACE;
    distance = numeric_limits<double>::max();
    final_position.resize(dimension_);
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        double current_distance;
        vector<double> current_position;
        
        if(surfaces_[i]->intersection(initial_position,
                                      initial_direction,
                                      current_distance,
                                      current_position)
           == Surface::Intersection::INTERSECTS)
        {
            
            if (current_distance < distance)
            {
                best_surface = i;
                distance = current_distance;
                final_position = current_position;
            }
        }
    }
    
    return best_surface;
}

// Only checks first intersection for each surface
int Solid_Geometry::
next_intersection(vector<double> const &initial_position,
                  vector<double> const &initial_direction,
                  int &final_region,
                  double &distance,
                  vector<double> &final_position) const
{
    int final_surface = NO_SURFACE;
    final_region = NO_REGION;
    distance = numeric_limits<double>::max();
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        double current_distance;
        vector<double> current_position;
        
        if(surfaces_[i]->intersection(initial_position,
                                      initial_direction,
                                      current_distance,
                                      current_position)
           == Surface::Intersection::INTERSECTS)
        {
            if (current_distance < distance)
            {
                vector<double> minus_position;
                vector<double> plus_position;
                
                new_position(-delta_distance(),
                             current_position,
                             initial_direction,
                             minus_position);
                new_position(delta_distance(),
                             current_position,
                             initial_direction,
                             plus_position);
                
                int minus_region = find_region(minus_position);
                int plus_region = find_region(plus_position);
                
                if (minus_region != plus_region)
                {
                    final_surface = i;
                    final_region = plus_region;
                    distance = current_distance;
                    final_position = current_position;
                }
            }
        }
    }
    
    return final_surface;
}

int Solid_Geometry::
next_boundary(vector<double> const &initial_position,
              vector<double> const &initial_direction,
              int &boundary_region,
              double &distance,
              vector<double> &final_position) const
{
    int surface = 0;
    double current_distance = numeric_limits<double>::max();
    vector<double> current_position(initial_position);
    distance = 0;

    while (surface != NO_SURFACE)
    {
        int temp_region = NO_REGION; // could be changed to plus_region
        
        surface = next_intersection(current_position,
                                    initial_direction,
                                    temp_region,
                                    current_distance,
                                    current_position);
        
        distance += current_distance;
        
        vector<double> minus_position;
        vector<double> plus_position;
        
        new_position(-delta_distance(),
                     current_position,
                     initial_direction,
                     minus_position);
        new_position(delta_distance(),
                     current_position,
                     initial_direction,
                     plus_position);
        
        int minus_region = find_region(minus_position);
        int plus_region = find_region(plus_position);
        
        if ((plus_region == NO_REGION || minus_region == NO_REGION) && (minus_region != plus_region))
        {
            boundary_region = (plus_region == NO_REGION) ? minus_region : plus_region;
            final_position = current_position;
            
            break;
        }
        
    }
    
    return surface;
}

void Solid_Geometry::
new_position(double distance,
             vector<double> const &initial_position,
             vector<double> const &initial_direction,
             vector<double> &final_position) const
{
    final_position.resize(dimension_);

    for (int i = 0; i < dimension_; ++i)
    {
        final_position[i] = initial_position[i] + initial_direction[i] * distance;
    }
}
