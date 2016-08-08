#ifndef Spatial_Discretization_Parser_hh
#define Spatial_Discretization_Parser_hh

#include "Parser.hh"
#include "Spatial_Discretization.hh"

class Finite_Element_Mesh;
class RBF_Mesh;
class Solid_Geometry;

/*
  Create a Spatial_Discretization object from XML file
*/
class Spatial_Discretization_Parser : public Parser<Spatial_Discretization>
{
public:
    
    // Creator
    Spatial_Discretization_Parser(pugi::xml_node &input_file);
    
    // Return pointer to object
    virtual shared_ptr<Spatial_Discretization> get_ptr() override
    {
        return spatial_;
    }
    
    // Parse a finite element mesh
    shared_ptr<Finite_Element_Mesh> get_fem(pugi::xml_node &spatial);
    
    // Parse a radial basis function mesh
    shared_ptr<RBF_Mesh> get_rbf_1d(pugi::xml_node &spatial);
    shared_ptr<RBF_Mesh> get_rbf_solid(pugi::xml_node &spatial);
    
private:
    
    shared_ptr<Spatial_Discretization> spatial_;
};

#endif
