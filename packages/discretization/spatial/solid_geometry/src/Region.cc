#include "Region.hh"

#include "Check.hh"

Region::
Region():
    material_(-1), // material will be >=0 once region is initialized
    max_recursion_(100)
{
}

void Region::
initialize(int material,
           vector<Surface::Relation> const &surface_relations,
           vector<shared_ptr<Surface> > const &surfaces,
           vector<Region::Relation> const &region_relations,
           vector<shared_ptr<Region> > const &regions)
{
    material_ = material;
    surface_relations_ = surface_relations;
    surfaces_ = surfaces;
    region_relations_ = region_relations;
    regions_ = regions;

    for (int i = 0; i < surfaces_.size(); ++i)
    {
        Assert(surfaces_[i]);
    }

    for (int i = 0; i < regions_.size(); ++i)
    {
        Assert(regions_[i]);
    }
}

Region::Relation Region::
relation(vector<double> const &point,
         int recursion_level) const
{
    if (recursion_level > max_recursion_)
    {
        AssertMsg(false, "recursion limit reached");
    }
    
    Surface::Relation surface_relation;
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        surface_relation = surfaces_[i]->relation(point);
        
        if (surface_relation != surface_relations_[i])
        {
            return Region::Relation::OUTSIDE;
        }
    }
    
    Region::Relation region_relation;
    
    for (int i = 0; i < regions_.size(); ++i)
    {
        region_relation = regions_[i]->relation(point,
                                                recursion_level + 1);
        
        if (region_relation != region_relations_[i])
        {
            return Region::Relation::OUTSIDE;
        }
    }
    
    return Region::Relation::INSIDE;
}
