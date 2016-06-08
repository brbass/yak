#include "Region.hh"

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
    regions_ = regions;
}

Region::Relation Region::
relation(vector<double> const &point,
         int recursion_level)
{
    if (recursion_level_ > max_recursion_)
    {
        AssertMsg(false, "recursion limit reached");
    }
    
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
        region_relation = regions_[i]->relation(point,
                                                recursion_level + 1);
        
        if (surface_relation != surface_relations_[i])
        {
            return Region::OUTSIDE;
        }
    }
    
    return Region::INSIDE;
}
