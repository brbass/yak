#ifndef Solver_Parser_hh
#define Solver_Parser_hh

#include "DFEM_Sweep_1D.hh"
#include "Discrete_To_Moment.hh"
#include "Fission.hh"
#include "Krylov_Iteration.hh"
#include "Local_RBF_Sweep_1D.hh"
#include "Moment_To_Discrete.hh"
#include "Ordinate_Sweep_Operator.hh"
#include "Parser.hh"
#include "RBF_Sweep_1D.hh"
#include "Scattering.hh"
#include "Solver.hh"
#include "Source_Iteration.hh"

/*
  Create a Solver object
*/
class Solver_Parser : public Parser<Solver>
{
public:

    // Constructor
    Solver_Parser(pugi::xml_node &input_file,
                  shared_ptr<Spatial_Discretization> spatial,
                  shared_ptr<Angular_Discretization> angular,
                  shared_ptr<Energy_Discretization> energy,
                  shared_ptr<Nuclear_Data> nuclear,
                  shared_ptr<Source_Data> source);

    // Return object
    virtual shared_ptr<Solver> get_ptr()
    {
        return solver_;
    }

    // Parse source iteration
    shared_ptr<Source_Iteration> parse_source_iteration();

    // Parse Krylov iteration
    shared_ptr<Krylov_Iteration> parse_krylov_iteration();

    // Parse the Ordinate_Sweep_Operator
    shared_ptr<Ordinate_Sweep_Operator> parse_sweeper();

    // Parse the DFEM sweeper
    shared_ptr<DFEM_Sweep_1D> parse_dfem();

    // Parse the RBF Sweeper
    shared_ptr<RBF_Sweep_1D> parse_rbf();

    // Parse the Local RBF Sweeper
    shared_ptr<Local_RBF_Sweep_1D> parse_rbf();
    
    // Parse the Discrete_To_Moment operator
    shared_ptr<Discrete_To_Moment> parse_discrete_to_moment();

    // Parse the Moment_To_Discrete operator
    shared_ptr<Moment_To_Discrete> parse_moment_to_discrete();

    // Parse the scattering operator
    shared_ptr<Scattering> parse_scattering();

    // Parse the fission operator
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
