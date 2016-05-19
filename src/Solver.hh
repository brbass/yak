#ifndef Solver_hh
#define Solver_hh

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

/*
  Pure virtual class for a solver of the transport equation
*/
class Solver
{
public:

    // Constructor
    Solver(shared_ptr<Spatial_Discretization> spatial_discretization,
           shared_ptr<Angular_Discretization> angular_discretization,
           shared_ptr<Energy_Discretization> energy_discretization,
           shared_ptr<Nuclear_Data> nuclear_data,
           shared_ptr<Source_Data> source_data);

    // Solve fixed source problem
    virtual void solve_steady_state(vector<double> &x) = 0;

    // Solve k-eigenvalue problem
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x) = 0;

    // Solve time-dependent problem
    virtual void solve_time_dependent(vector<double> &x) = 0;

    // Ouput data to XML file
    virtual void output(pugi::xml_node &output_node) const;

protected:

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    // shared_ptr<Temporal_Discretization> temporal_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
    shared_ptr<Source_Data> source_data_;
};

#endif
