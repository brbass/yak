#include "Solver_Parser.hh"

#include "Check.hh"
#include "DFEM_Sweep_1D.hh"
#include "Discrete_To_Moment.hh"
#include "Fission.hh"
#include "Krylov_Iteration.hh"
#include "Local_RBF_Sweep.hh"
#include "Moment_To_Discrete.hh"
#include "Null_Solver.hh"
#include "Sweep_Operator.hh"
#include "RBF_Sweep_1D.hh"
#include "Scattering.hh"
#include "Source_Iteration.hh"

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
    pugi::xml_node solver_node = input_file.child("solution_method");
    
    string solver_type = XML_Functions::child_value<string>(solver_node, "type");
    
    if (solver_type == "source_iteration")
    {
        solver_  = parse_source_iteration();
    }
    else if (solver_type == "krylov_iteration")
    {
        solver_ = parse_krylov_iteration();
    }
    else if (solver_type == "null_solver")
    {
        solver_ = make_shared<Null_Solver>(1,
                                           spatial,
                                           angular,
                                           energy,
                                           nuclear,
                                           source);
    }
    else
    {
        AssertMsg(false, "solver type " + solver_type + " not found");
    }
}    

shared_ptr<Source_Iteration> Solver_Parser::
parse_source_iteration()
{
    pugi::xml_node solver_node = input_file_.child("solution_method");
    
    int max_iterations = XML_Functions::child_value<int>(solver_node, "max_iterations");
    int solver_print = XML_Functions::child_value<int>(solver_node, "solver_print");
    double tolerance = XML_Functions::child_value<double>(solver_node, "tolerance");
    shared_ptr<Vector_Operator> sweeper = parse_sweeper();
    shared_ptr<Vector_Operator> discrete_to_moment = parse_discrete_to_moment();
    shared_ptr<Vector_Operator> moment_to_discrete = parse_moment_to_discrete();
    shared_ptr<Vector_Operator> scattering = parse_scattering();
    shared_ptr<Vector_Operator> fission = parse_fission();
    
    return make_shared<Source_Iteration>(max_iterations,
                                         solver_print,
                                         tolerance,
                                         spatial_,
                                         angular_,
                                         energy_,
                                         nuclear_,
                                         source_,
                                         sweeper,
                                         discrete_to_moment,
                                         moment_to_discrete,
                                         scattering,
                                         fission);
}

shared_ptr<Krylov_Iteration> Solver_Parser::
parse_krylov_iteration()
{
    pugi::xml_node solver_node = input_file_.child("solution_method");
    
    int max_iterations = XML_Functions::child_value<int>(solver_node, "max_iterations");
    int kspace = XML_Functions::child_value<int>(solver_node, "kspace");
    int solver_print = XML_Functions::child_value<int>(solver_node, "solver_print");
    double tolerance = XML_Functions::child_value<double>(solver_node, "tolerance");
    shared_ptr<Vector_Operator> sweeper = parse_sweeper();
    shared_ptr<Vector_Operator> discrete_to_moment = parse_discrete_to_moment();
    shared_ptr<Vector_Operator> moment_to_discrete = parse_moment_to_discrete();
    shared_ptr<Vector_Operator> scattering = parse_scattering();
    shared_ptr<Vector_Operator> fission = parse_fission();
    
    return make_shared<Krylov_Iteration>(max_iterations,
                                         kspace,
                                         solver_print,
                                         tolerance,
                                         spatial_,
                                         angular_,
                                         energy_,
                                         nuclear_,
                                         source_,
                                         sweeper,
                                         discrete_to_moment,
                                         moment_to_discrete,
                                         scattering,
                                         fission);
}

shared_ptr<Sweep_Operator> Solver_Parser::
parse_sweeper()
{
    pugi::xml_node solver = input_file_.child("solution_method");
    pugi::xml_node sweeper = solver.child("sweeper");
    
    string sweeper_type = XML_Functions::child_value<string>(sweeper, "type");
    
    if (sweeper_type == "dfem")
    {
        return parse_dfem();
    }
    else if (sweeper_type == "rbf")
    {
        return parse_rbf();
    }
    else if (sweeper_type == "rbf_local")
    {
        return parse_rbf_local();
    }
    else
    {
        AssertMsg(false, "sweeper type " + sweeper_type + " not found");
        
        return shared_ptr<Sweep_Operator>();
    }
}

shared_ptr<DFEM_Sweep_1D> Solver_Parser::
parse_dfem()
{
    return make_shared<DFEM_Sweep_1D>(spatial_,
                                      angular_,
                                      energy_,
                                      nuclear_,
                                      source_);
}

shared_ptr<RBF_Sweep_1D> Solver_Parser::
parse_rbf()
{
    return make_shared<RBF_Sweep_1D>(spatial_,
                                     angular_,
                                     energy_,
                                     nuclear_,
                                     source_);
}

shared_ptr<Local_RBF_Sweep> Solver_Parser::
parse_rbf_local()
{
    return make_shared<Local_RBF_Sweep>(spatial_,
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

