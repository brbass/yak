#ifndef Spatial_Discretization_Parser_hh
#define Spatial_Discretization_Parser_hh

#include "Spatial_Discretization.hh"

class Spatial_Discretization_Parser : public Parser<Spatial_Discretization>
{
public:
    
    Spatial_Discretization_Parser(pugi::xml_node &input_file);
    
    virtual shared_ptr<Spatial_Discretization> get_ptr()
    {
        return spatial_;
    }

    shared_ptr<Finite_Element_Mesh> get_fem(pugi::xml_node &spatial);
    shared_ptr<RBF_Mesh> get_rbf(pugi::xml_node &spatial);
    
private:
    
    shared_ptr<Spatial_Discretization> spatial_;
};

#endif
