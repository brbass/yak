#include "Solver_Parser.hh"

using namespace std;

Solver_Parser::
Solver_Parser(pugi::xml_node &input_file,
              shared_ptr<Spatial_Discretization> spatial,
              shared_ptr<Angular_Discretization> angular,
              shared_ptr<Energy_Discretization> energy,
              shared_ptr<Nuclear_Data> nuclear,
              shared_ptr<Source_Data> source):
    Parser(input_file),
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    nuclear_(nuclear),
    source_(source)
{
    pugi::xml_node solver_node = input_file.child("solver");
    
    string solver_type = child_value<string>(solver_node, "solver_type");
    
    if (solver_type == "source_iteration")
    {
        solver_  = parse_source_iteration();
    }
    else
    {
        AssertMsg(false, "solver type " + solver_type + " not found");
    }
}    

shared_ptr<Source_Iteration> Solver_Parser::
parse_source_iteration()
{
    pugi::xml_node solver_node = input_file.child("solver");
    
    int max_iterations = child_value<int>(solver_node, "max_iterations");
    double tolerance = child_value<int>(solver_node, "tolerance");
    shared_ptr<Vector_Operator> sweeper = parse_sweeper();
    shared_ptr<Vector_Operator> discrete_to_moment = parse_discrete_to_moment();
    shared_ptr<Vector_Operator> moment_to_discrete = parse_moment_to_discrete();
    shared_ptr<Vector_Operator> scattering = parse_scattering();
    shared_ptr<Vector_Operator> fission = parse_fission();
    
    return make_shared<Source_Iteration>(max_iterations,
                                         tolerance,
                                         spatial,
                                         angular,
                                         energy,
                                         nuclear,
                                         source,
                                         sweeper,
                                         discrete_to_moment,
                                         moment_to_discrete,
                                         scattering,
                                         fission);
}

shared_ptr<Ordinate_Sweep_Operator> Solver_Parser::
parse_sweeper()
{
    pugi::xml_node sweeper = input_file_.child("sweeper");
    
    string sweeper_type = child_value<string>(sweeper, "type");

    if (sweeper_type == "DFEM")
    {
        return parse_dfem();
    }
    else
    {
        InsistMsg(false, "sweeper type " + sweeper_type + " not found");
        
        return shared_ptr<Ordinate_Sweep_Operator>();
    }
}

shared_ptr<DFEM_Sweep> Solver_Parser::
parse_dfem()
{
    return make_shared<DFEM_Sweep>(spatial_,
                                   angular_,
                                   energy_,
                                   nuclear_,
                                   source_);
}

shared_ptr<Discrete_To_Moment> Solver_Parser::
parse_discrete_to_moment()
{
    return make_shared<Discrete_To_Moment>(spatial_,
                                           angular_,
                                           energy_);
}

shared_ptr<Moment_To_Discrete> Solver_Parser::
parse_moment_to_discrete()
{
    return make_shared<Moment_To_Discrete>(spatial_,
                                           angular_,
                                           energy_);
}

shared_ptr<Scattering> Solver_Parser::
parse_scattering()
{
    return make_shared<Scattering>(spatial_,
                                   angular_,
                                   energy_,
                                   nuclear_);
}

shared_ptr<Fission> Solver_Parser::
parse_fission()
{
    return make_shared<Fission>(spatial_,
                                angular_,
                                energy_,
                                nuclear_);
}
