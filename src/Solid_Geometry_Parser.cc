#include "Solid_Geometry_Parser.hh"

#include <iterator>

using namespace std;

Solid_Geometry_Parser::
Solid_Geometry_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node solid_node = input_file.child("solid_geometry");

    solid_ = get_solid(solid_node);
}

shared_ptr<Solid_Geometry> Solid_Geometry_Parser::
get_solid(pugi::xml_node &solid_node)
{
    int dimension = child_value<int>(solid_node, "dimension");
    
    pugi::xml_node surfaces_node = solid_node.child("surfaces");
    vector<shared_ptr<Surface> > surfaces = get_surfaces(surfaces_node,
                                                         dimension);

    pugi::xml_node regions_node = solid_node.child("regions");
    vector<shared_ptr<Region> > regions = get_regions(regions_node,
                                                      surfaces);
    
    return make_shared<Solid_Geometry>(dimension,
                                       surfaces,
                                       regions);
}

vector<shared_ptr<Surface> > Solid_Geometry_Parser::
get_surfaces(pugi::xml_node &surfaces_node,
             int dimension)
{
    int number_of_surfaces = distance(surfaces_node.children("surface").begin(),
                                      surfaces_node.children("surface").end());
    
    vector<shared_ptr<Surface> > surfaces(number_of_surfaces);

    for (pugi::xml_node surface_node = surfaces_node.child("surface"); surface_node; surface_node = surface_node.next_sibling("surface"))
    {
        int number = child_value<int>(surface_node, "number");
        string shape = child_value<string>(surface_node, "shape");
        
        shared_ptr<Surface> surface;

        Surface::Surface_Type surface_type;

        string type = child_value<int>(surface_node, "type");

        if (type == "vacuum_boundary")
        {
            surface_type = Surface::Surface_Type::VACUUM_BOUNDARY;
        }
        else if (type == "reflective_boundary")
        {
            surface_type = Surface::Surface_Type::REFLECTIVE_BOUNDARY;
        }
        else if (type == "internal")
        {
            surface_type = Surface::Surface_Type::INTERNAL;
        }
        else
        {
            AssertMsg(false, "surface type not found");
        }
        
        if (shape == "line")
        {
            vector<double> origin = child_vector<double>(surface_node, "origin");
            vector<double> normal = child_vector<double>(surface_node, "normal");
            
            surface = make_shared<Line>(surface_type,
                                        origin,
                                        normal);
        }
        else if (shape == "circle")
        {
            double radius = child_value<double>(surface_node, "radius");
            vector<double> origin = child_vector<double>(surface_node, "origin");

            surface = make_shared<Circle>(surface_type,
                                          radius,
                                          origin);
        }
        else if (shape == "plane")
        {
            vector<double> origin = child_vector<double>(surface_node, "origin");
            vector<double> normal = child_vector<double>(surface_node, "normal");
            
            surface = make_shared<Line>(surface_type,
                                        origin,
                                        normal);
        }
        else if (shape == "sphere")
        {
            double radius = child_value<double>(surface_node, "radius");
            vector<double> origin = child_vector<double>(surface_node, "origin");

            surface = make_shared<Sphere>(surface_type,
                                          radius,
                                          origin);
        }
        else if (shape == "cylinder")
        {
            double radius = child_value<double>(surface_node, "radius");
            vector<double> origin = child_vector<double>(surface_node, "origin");
            vector<double> direction = child_vector<double>(surface_node, "direction");
            
            surface = make_shared<Sphere>(surface_type,
                                          radius,
                                          origin,
                                          direction);
            
        }
        else
        {
            AssertMsg(false, "surface shape not found");
        }
        
        surfaces[number] = surface;
    }
    
    return surfaces;
}

vector<shared_ptr<Region> > Solid_Geometry_Parser::
get_regions(pugi::xml_node &regions_node,
            vector<shared_ptr<Surface> > const &surfaces)
{
    int number_of_regions = distance(regions_node.children("region").begin(),
                                     regions_node.children("region").end());
    
    vector<shared_ptr<Region> > regions(number_of_regions);

    for (int i = 0; i < number_of_regions; ++i)
    {
        regions[i] = make_shared<Region>();
    }
    
    for (pugi::xml_node region_node = regions_node.child("region"); region_node; region_node = region_node.next_sibling("region"))
    {
        int number = child_value<int>(region_node, "number");
        int material = child_value<int>(region_node, "material");

        vector<Surface::Relation> surface_relations;
        vector<shared_ptr<Surface> > local_surfaces;
 
        for (pugi::xml_node surface_relation_node = region_node.child("surface_relation"); surface_relation_node; surface_relation_node = surface_relation_node.next_sibling("surface_relation"))
        {
            int relation_number = child_value<int>(surface_relation_node, "surface");
            string relation_string = child_value<string>(surface_relation_node, "relation");
            
            Surface::Relation relation;

            if (relation_string == "positive")
            {
                relation = Surface::Relation::POSITIVE;
            }
            else if (relation_string == "negative")
            {
                relation = Surface::Relation::NEGATIVE;
            }
            else if (relation_string == "equal")
            {
                relation = Surface::Relation::EQUAL;
            }
            else if (relation_string == "outside")
            {
                relation = Surface::Relation::OUTSIDE;
            }
            else if (relation_string == "inside")
            {
                relation = Surface::Relation::INSIDE;
            }
            else
            {
                AssertMsg(false, "surface relation not found");
            }

            surface_relations.push_back(relation);
            local_surfaces.push_back(surfaces[relation_number]);
        }       
        
        vector<Region::Relation> region_relations;
        vector<shared_ptr<Region> > local_regions;
        
        for (pugi::xml_node region_relation_node = region_node.child("region_relation"); region_relation_node; region_relation_node = region_relation_node.next_sibling("region_relation"))
        {
            int relation_number = child_value<int>(region_relation_node, "region");
            string relation_string = child_value<string>(region_relation_node, "relation");
            
            Region::Relation relation;

            if (relation_string == "outside")
            {
                relation = Region::Relation::OUTSIDE;
            }
            else if (relation_string == "inside")
            {
                relation = Region::Relation::INSIDE;
            }
            else
            {
                AssertMsg(false, "region relation not found");
            }

            region_relations.push_back(relation);
            local_regions.push_back(regions[relation_number]);
        }       
        
        regions[number]->initialize(material,
                                    surface_relations,
                                    local_surfaces,
                                    region_relations,
                                    local_regions);
        
    }
    
    return regions;
}
