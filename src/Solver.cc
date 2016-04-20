#include "Solver.hh"

Solver::
Solver(shared_ptr<Spatial_Discretization> spatial_discretization,
       shared_ptr<Angular_Discretization> angular_discretization,
       shared_ptr<Energy_Discretization> energy_discretization,
       shared_ptr<Nuclear_Data> nuclear_data,
       shared_ptr<Source_Data> source_data):
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    nuclear_data_(nuclear_data),
    source_data_(source_data)
{
}
