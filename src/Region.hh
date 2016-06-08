#ifndef Region_hh
#define Region_hh

#include "Surface.hh"

/* 
   Describes a solid geometry region
   
   The memory for the regions must be allocated before the data is initialized,
   as a region can be defined in terms of other regions. Because of this, the
   "initialize" function must be called after the constructor. This allows 
   construction of the regions in any order without regard for interdependencies.
*/
class Region
{
public:

    // Is the point inside or outside of the region
    enum class Relation
    {
        INSIDE,
        OUTSIDE
    };

    // Constructor
    Region();

    // Initialization values
    void initialize(int material,
                    vector<Surface::Relation> const &surface_relations,
                    vector<shared_ptr<Surface> > const &surfaces,
                    vector<Region::Relation> const &region_relations,
                    vector<shared_ptr<Region> > const &regions);
    
    int material() const
    {
        return material_;
    }
    Surface::Relation surface_relation(int s) const
    {
        return surface_relations_[s];
    }
    shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    Relation region_relation(int r) const
    {
        return region_relations_[r];
    }
    shared_ptr<Region> const &region(int r) const
    {
        return region_[r];
    }
    
    Relation relation(vector<double> const &point,
                      int recursion_level = 0) const;
    
private:
    
    int material_;
    int max_recursion_;
    vector<Surface::Relation> surface_relations_;
    vector<shared_ptr<Surface> > surfaces_;
    vector<Region::Relation> region_relations_;
    vector<shared_ptr<Region> > regions_;
}

#endif
