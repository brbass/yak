#ifndef Spatial_Discretization_Parser_hh
#define Spatial_Discretization_Parser_hh

#include "Finite_Element_Mesh.hh"
#include "Local_RBF_Mesh.hh"
#include "Parser.hh"
#include "RBF_Mesh.hh"
#include "Spatial_Discretization.hh"

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

    void get_solid_points(shared_ptr<Solid_Geometry> solid_geometry,
                          pugi::xml_node &spatial,
                          int &number_of_points,
                          int &number_of_boundary_points,
                          int &number_of_internal_points,
                          vector<int> &material,
                          vector<int> &boundary_points,
                          vector<int> &internal_points,
                          vector<double> &positions,
                          vector<double> &boundary_normal);
    
    shared_ptr<Spatial_Discretization> spatial_;
};

#endif
