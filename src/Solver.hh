#ifndef Solver_hh
#define Solver_hh

#include <memory>
#include <string>
#include <vector>

#include "pugixml.hpp"

class Angular_Discretization;
class Energy_Discretization;
class Nuclear_Data;
class Source_Data;
class Spatial_Discretization;

using std::shared_ptr;
using std::string;
using std::vector;

/*
  Pure virtual class for a solver of the transport equation
*/
class Solver
{
public:

    // Constructor
    Solver(int solver_print,
           shared_ptr<Spatial_Discretization> spatial_discretization,
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

    virtual void print_name(string solution_type);
    virtual void print_iteration(int iteration);
    virtual void print_convergence();
    virtual void print_failure();

    int solver_print_;
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    // shared_ptr<Temporal_Discretization> temporal_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
    shared_ptr<Source_Data> source_data_;
};

#endif
