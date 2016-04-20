#include "Driver.hh"

#include "pugixml.hpp"

#include "Angular_Discretization_Parser.hh"
#include "Check.hh"
#include "Energy_Discretization_Parser.hh"
#include "Nuclear_Data_Parser.hh"
#include "Solver_Parser.hh"
#include "Source_Data_Parser.hh"
#include "Sweeper_Parser.hh"
#include "Spatial_Discretization_Parser.hh"
#include "Transport_Problem_Parser.hh"
#include "Vector_Operator.hh"
#include "XML_Child_Value.hh"

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
    
    if (!input_file.load_file(xml_in_.c_str()))
    {
        AssertMsg(false, "Could not open xml input file \"" + xml_in_ + "\"");
    }
    
    pugi::xml_node input_file = input_document.child("input");
    
    // level 1 classes

    Spatial_Discretization_Parser spatial_parser(input_file_);
    Angular_Discretization_Parser angular_parser(input_file_);
    Energy_Discretization_Parser energy_parser(input_file_);
    // Temporal_Discretization_Parser temporal_parser(input_file_);
    
    shared_ptr<Spatial_Discretization> spatial = spatial_parser.get_ptr();
    shared_ptr<Angular_Discretization> angular = angular_parser.get_ptr();
    shared_ptr<Energy_Discretization> energy = energy_parser.get_ptr();
    // shared_ptr<Temporal_Discretization> temporal = temporal_parser.get_ptr();
    
    // level 2 classes
    
    Nuclear_Data_Parser nuclear_parser(input_file_,
                                       spatial,
                                       angular,
                                       energy);
    Source_Data_Parser source_parser(input_file_,
                                     spatial,
                                     angular,
                                     energy);
    
    shared_ptr<Nuclear_Data> nuclear = nuclear_parser.get_ptr();
    shared_ptr<Source_Data> source = source_parser.get_ptr();
    
    // level 3 class
    
    Solver_Parser solver_parser(input_file_,
                                spatial,
                                angular,
                                energy,
                                nuclear,
                                source);
    
    shared_ptr<Solver> solver = solver_parser.get_ptr();
    
    // level 4 class
    
    string transport_type = child_value<string>(input_file_, "transport_type");
    
    Transport_Problem_Parser transport_parser(input_file_,
                                              solver);
    
    shared_ptr<Transport_Problem> transport = transport_parser.get_ptr();

    // solve problem
    
    problem.solve();
    
    pugi::xml_document output_document;
    pugi::xml_node output_file = output_document.append_child("output");
    
    problem.output(xml_out_);
}
