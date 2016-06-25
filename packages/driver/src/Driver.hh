#ifndef Driver_hh
#define Driver_hh

#include <string>

using std::string;

/*
  Create and run transport problem from XML file
*/
class Driver
{
public:

    // Constructor
    Driver(string filename);

private:

    // Run transport problem
    void run_problem();
    
    string xml_in_;
    string xml_out_;
};

#endif
