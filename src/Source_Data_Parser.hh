#ifndef Source_Data_Parser_hh
#define Source_Data_Parser_hh

#include "Parser.hh"
#include "Source_Data.hh"

class Source_Data_Parser : public Parser<Source_Data>
{
public:
    
    Source_Data_Parser(pugi::xml_node &input_file,
                       shared_ptr<Spatial_Discretization> spatial,
                       shared_ptr<Angular_Discretization> angular,
                       shared_ptr<Energy_Discretization> energy);
    
    virtual shared_ptr<Source_Data> get_ptr()
    {
        return source_;
    }

private:
    
    shared_ptr<Source_Data> source_;
    
    shared_ptr<Spatial_Discretization> spatial_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
};

#endif
