#ifndef Null_Solver_hh
#define Null_Solver_hh

#include "Solver.hh"

/*
  Allows setup but not solution of problem
*/
class Null_Solver : public Solver
{
public:

    Null_Solver(int solver_print,
                shared_ptr<Spatial_Discretization> spatial_discretization,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                shared_ptr<Nuclear_Data> nuclear_data,
                shared_ptr<Source_Data> source_data);

    virtual void solve_steady_state(vector<double> &x) override {};
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x) override {};
    virtual void solve_time_dependent(vector<double> &x) override {};
    virtual void output(pugi::xml_node &output_node) const override;
};

#endif
