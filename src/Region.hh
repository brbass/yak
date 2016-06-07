#ifndef Region_hh
#define Region_hh

#include "Surface.hh"

class Region
{
public:

    enum class Relation
    {
        INSIDE,
        OUTSIDE
    };
    
    Region(int material,
           vector<Surface::Relation> surface_relations,
           vector<shared_ptr<Surface> > surfaces,
           vector<Region::Relation> region_relations,
           vector<shared_ptr<Region> > regions);
    
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
    
    Relation relation(vector<double> &point) const;
    
private:

    int material_;
    vector<Surface::Relation> surface_relations_;
    vector<shared_ptr<Surface> > surfaces_;
    vector<Region::Relation> region_relations_;
    vector<shared_ptr<Region> > regions_;
}

#endif
