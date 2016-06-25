#include "Null_Solver.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

Null_Solver::
Null_Solver(int solver_print,
            shared_ptr<Spatial_Discretization> spatial_discretization,
            shared_ptr<Angular_Discretization> angular_discretization,
            shared_ptr<Energy_Discretization> energy_discretization,
            shared_ptr<Nuclear_Data> nuclear_data,
            shared_ptr<Source_Data> source_data):
    Solver(solver_print,
           spatial_discretization,
           angular_discretization,
           energy_discretization,
           nuclear_data,
           source_data)
{
}
    
void Null_Solver::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("krylov_iteration");
    
    nuclear_data_->output(output_node);
    source_data_->output(output_node);
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}
