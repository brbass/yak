#ifndef Solver_Parser_hh
#define Solver_Parser_hh

#include "DFEM_Sweep.hh"
#include "Discrete_To_Moment.hh"
#include "Fission.hh"
#include "Moment_To_Discrete.hh"
#include "Ordinate_Sweep_Operator.hh"
#include "Parser.hh"
#include "Scattering.hh"
#include "Solver.hh"
#include "Source_Iteration.hh"

class Solver_Parser : public Parser<Solver>
{
public:
    
    Solver_Parser(pugi::xml_node &input_file,
                  shared_ptr<Spatial_Discretization> spatial,
                  shared_ptr<Angular_Discretization> angular,
                  shared_ptr<Energy_Discretization> energy,
                  shared_ptr<Nuclear_Data> nuclear,
                  shared_ptr<Source_Data> source);
    
    virtual shared_ptr<Solver> get_ptr()
    {
        return solver_;
    }
    
    shared_ptr<Source_Iteration> parse_source_iteration();
    shared_ptr<Ordinate_Sweep_Operator> parse_sweeper();
    shared_ptr<DFEM_Sweep> parse_dfem();
    shared_ptr<Discrete_To_Moment> parse_discrete_to_moment();
    shared_ptr<Moment_To_Discrete> parse_moment_to_discrete();
    shared_ptr<Scattering> parse_scattering();
    shared_ptr<Fission> parse_fission();

private:
    
    shared_ptr<Solver> solver_;

    shared_ptr<Spatial_Discretization> spatial_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
    shared_ptr<Nuclear_Data> nuclear_;
    shared_ptr<Source_Data> source_;
};

#endif
