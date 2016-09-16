#include "Driver.hh"

#include "pugixml.hh"

#include "Angular_Discretization_Parser.hh"
#include "Check.hh"
#include "Energy_Discretization_Parser.hh"
#include "Nuclear_Data_Parser.hh"
#include "Solver_Parser.hh"
#include "Source_Data_Parser.hh"
#include "Spatial_Discretization_Parser.hh"
#include "Timer.hh"
#include "Transport_Problem_Parser.hh"
#include "Vector_Operator.hh"
#include "XML_Functions.hh"

Driver::
Driver(string xml_in):
    xml_in_(xml_in)
{
    xml_out_ = xml_in_ + ".out";
    
    run_problem();
}

void Driver::
run_problem()
{
    pugi::xml_document input_document;
    
    if (!input_document.load_file(xml_in_.c_str()))
    {
        AssertMsg(false, "Could not open xml input file \"" + xml_in_ + "\"");
    }
    
    pugi::xml_node input_file = input_document.child("input");
    
    Timer total_timer;
    total_timer.start();

    Timer timer;
    
    // level 1 classes

    timer.start();
        
    Spatial_Discretization_Parser spatial_parser(input_file);
    Angular_Discretization_Parser angular_parser(input_file);
    Energy_Discretization_Parser energy_parser(input_file);
    // Temporal_Discretization_Parser temporal_parser(input_file);
    
    shared_ptr<Spatial_Discretization> spatial = spatial_parser.get_ptr();
    shared_ptr<Angular_Discretization> angular = angular_parser.get_ptr();
    shared_ptr<Energy_Discretization> energy = energy_parser.get_ptr();
    // shared_ptr<Temporal_Discretization> temporal = temporal_parser.get_ptr();

    timer.stop();
    discretization_parser_time_ = timer.time();
    
    // level 2 classes

    timer.start();
    
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

    timer.stop();
    data_parser_time_ = timer.time();
    
    // level 3 class

    timer.start();
    
    Solver_Parser solver_parser(input_file,
                                spatial,
                                angular,
                                energy,
                                nuclear,
                                source);
    
    shared_ptr<Solver> solver = solver_parser.get_ptr();

    timer.stop();
    solver_parser_time_ = timer.time();
    
    // level 4 class

    timer.start();
    
    Transport_Problem_Parser transport_parser(input_file,
                                              solver);
    
    shared_ptr<Transport_Problem> transport = transport_parser.get_ptr();

    timer.stop();
    transport_parser_time_ = timer.time();
    
    // solve problem

    timer.start();
    
    transport->solve();

    timer.stop();
    solution_time_ = timer.time();

    // output data

    pugi::xml_document output_document;
    pugi::xml_node output_file = output_document.append_child("output");
    
    transport->output(output_file);
    
    total_timer.stop();
    total_time_ = total_timer.time();
    output(output_file);
    
    output_document.save_file(xml_out_.c_str());
}

void Driver::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node timing = output_node.append_child("timing");
    
    XML_Functions::append_child(timing, total_time_, "total");
    XML_Functions::append_child(timing, discretization_parser_time_, "discretization_parser");
    XML_Functions::append_child(timing, data_parser_time_, "data_parser");
    XML_Functions::append_child(timing, solver_parser_time_, "solver_parser");
    XML_Functions::append_child(timing, transport_parser_time_, "transport_parser");
    XML_Functions::append_child(timing, solution_time_, "solution");
}
