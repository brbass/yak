#ifndef Angular_Discretization_Parser_hh
#define Angular_Discretization_Parser_hh

#include "Angular_Discretization.hh"

class Angular_Discretization_Parser : public Parser<Angular_Discretization>
{
public:
    
    Angular_Discretization_Parser(pugi::xml_node &input_file);
    
    virtual shared_ptr<Angular_Discretization> get_ptr()
    {
        return angular_;
    }
    
private:
    
    shared_ptr<Angular_Discretization> angular_;
};

#endif
