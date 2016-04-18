#ifndef Transport_Problem_Parser_hh
#define Transport_Problem_Parser_hh

#include "Transport_Problem.hh"

class Transport_Problem_Parser
{
public:
    
    Transport_Problem_Parser(pugi::xml_node input_file);

    unique_ptr<Transport_Problem> parse();

private:
    
};

#endif
