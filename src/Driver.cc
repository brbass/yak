#include "Driver.hh"

#include "pugixml.hpp"

#include "Angular_Discretization_Parser.hh"
#include "Check.hh"
#include "Energy_Discretization_Parser.hh"
#include "Nuclear_Data_Parser.hh"
#include "Solver_Parser.hh"
#include "Source_Data_Parser.hh"
#include "Spatial_Discretization_Parser.hh"
#include "Transport_Problem_Parser.hh"
#include "Vector_Operator.hh"
// #include "XML_Functions.hh"

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
    // folder_ = filename.substr(0, filename.find_last_of("/\\") + 1);
    
    pugi::xml_document input_document;
    
    if (!input_document.load_file(xml_in_.c_str()))
    {
        AssertMsg(false, "Could not open xml input file \"" + xml_in_ + "\"");
    }
    
    pugi::xml_node input_file = input_document.child("input");
    
    // level 1 classes

    Spatial_Discretization_Parser spatial_parser(input_file);
    Angular_Discretization_Parser angular_parser(input_file);
    Energy_Discretization_Parser energy_parser(input_file);
    // Temporal_Discretization_Parser temporal_parser(input_file);
    
    shared_ptr<Spatial_Discretization> spatial = spatial_parser.get_ptr();
    shared_ptr<Angular_Discretization> angular = angular_parser.get_ptr();
    shared_ptr<Energy_Discretization> energy = energy_parser.get_ptr();
    // shared_ptr<Temporal_Discretization> temporal = temporal_parser.get_ptr();
    
    // level 2 classes
    
    Nuclear_Data_Parser nuclear_parser(input_file,
                                       spatial,
                                       angular,
                                       energy);
    Source_Data_Parser source_parser(input_file,
                                     spatial,
                                     angular,
                                     energy);
    
    shared_ptr<Nuclear_Data> nuclear = nuclear_parser.get_ptr();
    shared_ptr<Source_Data> source = source_parser.get_ptr();
    
    // level 3 class
    
    Solver_Parser solver_parser(input_file,
                                spatial,
                                angular,
                                energy,
                                nuclear,
                                source);
    
    shared_ptr<Solver> solver = solver_parser.get_ptr();
    
    // level 4 class
    
    Transport_Problem_Parser transport_parser(input_file,
                                              solver);
    
    shared_ptr<Transport_Problem> transport = transport_parser.get_ptr();

    // solve problem
    
    transport->solve();
    
    pugi::xml_document output_document;
    pugi::xml_node output_file = output_document.append_child("output");
    
    transport->output(output_file);
    
    output_document.save_file(xml_out_.c_str());
}
