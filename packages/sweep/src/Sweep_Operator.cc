#include "Sweep_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

namespace // anonymous
{
    int get_size(Sweep_Operator::Sweep_Type sweep_type,
                 shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Source_Data> source_data)
    {
        switch(sweep_type)
        {
        case Sweep_Operator::Sweep_Type::MOMENT:
            return (spatial_discretization->number_of_cells()
                    * spatial_discretization->number_of_nodes()
                    * energy_discretization->number_of_groups()
                    * angular_discretization->number_of_moments()
                    + source_data->number_of_augments());
        default: // Sweep_Type::ORDINATE
            return (spatial_discretization->number_of_cells()
                    * spatial_discretization->number_of_nodes()
                    * energy_discretization->number_of_groups()
                    * angular_discretization->number_of_ordinates()
                    + source_data->number_of_augments());
        }
    }
}

Sweep_Operator::
Sweep_Operator(Sweep_Type sweep_type,
               shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Nuclear_Data> nuclear_data,
               shared_ptr<Source_Data> source_data):
    Vector_Operator(get_size(sweep_type,
                             spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             source_data),
                    get_size(sweep_type,
                             spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             source_data)),
    include_boundary_source_(false),
    sweep_type_(sweep_type),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    nuclear_data_(nuclear_data),
    source_data_(source_data)
{
}

