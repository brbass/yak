#include "Scattering_Operator.hh"

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Scattering_Operator::
Scattering_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                    shared_ptr<Angular_Discretization> angular_discretization,
                    shared_ptr<Energy_Discretization> energy_discretization,
                    shared_ptr<Nuclear_Data> nuclear_data,
                    Scattering_Type scattering_type):
    Vector_Operator(get_size(spatial_discretization,
                             angular_discretization,
                             energy_discretization),
                    get_size(spatial_discretization,
                             angular_discretization,
                             energy_discretization)),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    nuclear_data_(nuclear_data),
    scattering_type_(scattering_type)
{
}

void Scattering_Operator::
apply(vector<double> &x)
{
    switch(scattering_type_)
    {
    case Scattering_Type::FULL:
        apply_full(x);
        break;
    case Scattering_Type::COHERENT:
        apply_coherent(x);
        break;
    case Scattering_Type::INCOHERENT:
        apply_incoherent(x);
        break;
    }
}

void Scattering_Operator::
apply_incoherent(vector<double> &x)
{
    vector<double> y(x);
    
    apply_full(y);
    
    apply_coherent(x);

    for (int i = 0; i < row_size(); ++i)
    {
        x[i] = y[i] - x[i];
    }
}

int Scattering_Operator::
get_size(shared_ptr<Spatial_Discretization> spatial_discretization,
         shared_ptr<Angular_Discretization> angular_discretization,
         shared_ptr<Energy_Discretization> energy_discretization)
{
    return (spatial_discretization->number_of_cells()
            * spatial_discretization->number_of_nodes()
            * energy_discretization->number_of_groups()
            * angular_discretization->number_of_moments());
}

