#include "Region.hh"

Region::
Region(int material,
       vector<Surface::Relation> surface_relations,
       vector<shared_ptr<Surface> > surfaces,
       vector<Region::Relation> region_relations,
       vector<shared_ptr<Region> > regions):
    material_(material),
    surface_relations_(surface_relations),
    surfaces_(surfaces),
    region_relations_(region_relations),
    regions_(regions)
{
}

Region::Relation Region::
relation(vector<double> &point)
{
    Surface::Relation surface_relation;
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        surface_relation = surfaces_[i]->relation(point);
        
        if (surface_relation != surface_relations_[i])
        {
            return Region::OUTSIDE;
        }
    }

    Region::Relation region_relation;
    
    for (int i = 0; i < regions_.size(); ++i)
    {
        region_relation = regions_[i]->relation(point);
        
        if (surface_relation != surface_relations_[i])
        {
            return Region::OUTSIDE;
        }
    }

    return Region::INSIDE;
}
