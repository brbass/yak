#include "Solver.hh"

#include <iostream>

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
print_name(string solution_type)
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
print_iteration(int iteration)
{
    if (solver_print_)
    {
        cout << "\t\titer:\t";
        cout << iteration;
        cout << endl;
    }
}

void Solver::
print_convergence()
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tConverged";
        cout << endl;
    }
}

void Solver::
print_failure()
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tFailed to converge";
        cout << endl;
    }
}
