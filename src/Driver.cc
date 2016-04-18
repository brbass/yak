#include "Driver.hh"

#include "pugixml.hpp"

#include "Check.hh"
#include "Transport_Problem_Parser.hh"

Driver::
Driver(string xml_in):
    xml_in_(xml_in)
{
    // xml_out_ = xml_in_.substr(0, filename.find_last_of(".")) + "-out.xml";
    xml_out_ = xml_in_ + ".out";
    
    run_problem();
}

void Driver::
run_problem()
{
    folder_ = filename.substr(0, filename.find_last_of("/\\") + 1);
    
    pugi::xml_document input_document;
    
    if (!input_file.load_file(xml_in_.c_str()))
    {
        AssertMsg(false, "Could not open input file \"" + xml_in_ + "\"");
    }
    
    pugi::xml_node input_file = input_document.child("input");
    
    Transport_Problem_Parser parser(input_file);
    
    unique_ptr<Transport_Problem> problem = parser.parse();
    
    problem.solve();
}
