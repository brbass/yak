#ifndef Solid_Geometry_Parser_hh
#define Solid_Geometry_Parser_hh

#include "Parser.hh"
#include "Solid_Geometry.hh"

class Solid_Geometry_Parser : public Parser<Solid_Geometry>
{
public:

    // Creator
    Solid_Geometry_Parser(pugi::xml_node &input_file);

    // Return pointer to object
    virtual shared_ptr<Solid_Geometry> get_ptr() override
    {
        return solid_;
    }

    
private:

    // Get solid geometry
    shared_ptr<Solid_Geometry> get_solid(pugi::xml_node &solid_node);
    vector<shared_ptr<Surface> > get_surfaces(pugi::xml_node &surfaces_node,
                                              int dimension);
    vector<shared_ptr<Region> > get_regions(pugi::xml_node &regions_node,
                                            vector<shared_ptr<Surface> > const &surfaces);
    
    shared_ptr<Solid_Geometry> solid_;
};

#endif
