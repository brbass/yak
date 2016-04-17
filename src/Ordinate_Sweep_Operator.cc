#include "Sweep_Operator.cc"

Sweep_Operator::
Sweep_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Nuclear_Data> nuclear_data,
               Scattering_Type scattering_type):
    Vector_Operator(spatial_discretization_->number_of_cells() * 
                    spatial_discretization_->number_of_nodes() * 
                    energy_discretization_->number_of_groups() * 
                    angular_discretization_->number_of_ordinates(),
                    spatial_discretization_->number_of_cells() * 
                    spatial_discretization_->number_of_nodes() * 
                    energy_discretization_->number_of_groups() * 
                    angular_discretization_->number_of_ordinates()),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    nuclear_data_(nuclear_data)
{
}
