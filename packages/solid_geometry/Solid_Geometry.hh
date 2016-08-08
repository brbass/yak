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

    enum Geometry_Errors
    {
        NO_REGION = -1,
        NO_SURFACE = -2
    };
    
    Solid_Geometry(int dimension,
                   vector<shared_ptr<Surface> > const &surfaces,
                   vector<shared_ptr<Region> > const &regions);
    
    int dimension() const
    {
        return dimension_;
    }
    int number_of_surfaces() const
    {
        return surfaces_.size();
    }
    int number_of_regions() const
    {
        return regions_.size();
    }
    int number_of_materials() const
    {
        return number_of_materials_;
    }
    double delta_distance() const
    {
        return delta_distance_;
    }
    
    shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    shared_ptr<Region> const &region(int r) const
    {
        return regions_[r];
    }
    
    int find_region(vector<double> const &position) const;
    
    int find_surface(vector<double> const &position) const;
    
    int next_intersection(vector<double> const &initial_position,
                          vector<double> const &initial_direction,
                          int &final_region,
                          double &distance,
                          vector<double> &final_position) const;

    int next_boundary(vector<double> const &initial_position,
                      vector<double> const &initial_direction,
                      int &boundary_region,
                      double &distance,
                      vector<double> &final_position) const;
    
    int next_geometric_intersection(vector<double> const &initial_position,
                                    vector<double> const &initial_direction,
                                    double &distance,
                                    vector<double> &final_position) const;
    
    void new_position(double distance,
                      vector<double> const &initial_position,
                      vector<double> const &initial_direction,
                      vector<double> &final_position) const;
    
private:

    int dimension_;
    int number_of_materials_;
    double delta_distance_;
    vector<shared_ptr<Surface> > surfaces_;
    vector<shared_ptr<Region> > regions_;
};

#endif
