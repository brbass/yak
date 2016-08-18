#ifndef Driver_hh
#define Driver_hh

#include <string>

#include "pugixml.hh"

using std::string;

/*
  Create and run transport problem from XML file
*/
class Driver
{
public:

    // Constructor
    Driver(string filename);

    void output(pugi::xml_node &output_node) const;
    
private:

    // Run transport problem
    void run_problem();

    double total_time_;
    double discretization_parser_time_;
    double data_parser_time_;
    double solver_parser_time_;
    double transport_parser_time_;
    double solution_time_;
    
    string xml_in_;
    string xml_out_;
};

#endif
