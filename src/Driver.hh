#include <string>

using std::string;

class Driver
{
public:
    
    Driver(string filename);

private:
    
    void run_problem();
    
    string xml_in_;
    string xml_out_;
};
