#ifndef Nuclear_Data_Parser_hh
#define Nuclear_Data_Parser_hh

#include "Nuclear_Data.hh"
#include "Parser.hh"

/*
  Create a Nuclear_Data object from XML file
*/
class Nuclear_Data_Parser : public Parser<Nuclear_Data>
{
public:

    // Constructor
    Nuclear_Data_Parser(pugi::xml_node &input_file,
                        shared_ptr<Spatial_Discretization> spatial,
                        shared_ptr<Angular_Discretization> angular,
                        shared_ptr<Energy_Discretization> energy);

    // Return Nuclear_Data object
    virtual shared_ptr<Nuclear_Data> get_ptr()
    {
        return nuclear_;
    }

private:
    
    shared_ptr<Nuclear_Data> nuclear_;

    shared_ptr<Spatial_Discretization> spatial_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
};

#endif
