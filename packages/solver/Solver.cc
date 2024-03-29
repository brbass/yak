#include "Solver.hh"

#include <iomanip>
#include <iostream>
#include <string>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Solver::
Solver(int solver_print,
       shared_ptr<Spatial_Discretization> spatial_discretization,
       shared_ptr<Angular_Discretization> angular_discretization,
       shared_ptr<Energy_Discretization> energy_discretization,
       shared_ptr<Nuclear_Data> nuclear_data,
       shared_ptr<Source_Data> source_data):
    solver_print_(solver_print),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    nuclear_data_(nuclear_data),
    source_data_(source_data)
{
}

void Solver::
output(pugi::xml_node &output_node) const
{
    nuclear_data_->output(output_node);
    source_data_->output(output_node);
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}

void Solver::
print_name(string solution_type) const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\t*******************************************************";
        cout << endl;
        cout << "\t\t***** " << solution_type;
        cout << endl;
        cout << "\t\t*******************************************************";
        cout << endl;
    }        
}

void Solver::
print_iteration(int iteration) const
{
    if (solver_print_)
    {
        cout << "\t\titer:\t";
        cout << iteration;
        cout << "\t";
    }
}

void Solver::
print_convergence() const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tConverged";
        cout << endl;
    }
}

void Solver::
print_failure() const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tFailed to converge";
        cout << endl;
    }
}

void Solver::
print_value(double value) const
{
    if (solver_print_)
    {
        cout << "value:\t";
        cout << value;
        cout << "\t";
    }
}

void Solver::
print_error(double error) const
{
    if (solver_print_)
    {
        cout << "error:\t";
        cout << error;
        cout << endl;
    }
}

void Solver::
print_eigenvalue(double eigenvalue) const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tk_eigenvalue:\t";
        cout << setprecision(10);
        cout << eigenvalue;
        cout << endl;
    }
}
