#ifndef Energy_Discretization_Parser_hh
#define Energy_Discretization_Parser_hh

#include "Energy_Discretization.hh"

class Energy_Discretization_Parser : public Parser<Energy_Discretization>
{
public:
    
    Energy_Discretization_Parser(pugi::xml_node &input_file);
    
    virtual shared_ptr<Energy_Discretization> get_ptr()
    {
        return energy_;
    }
    
private:
    
    shared_ptr<Energy_Discretization> energy_;
};

#endif
