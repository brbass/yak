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
    Vector_Operator(spatial_discretization->number_of_cells() * 
                    spatial_discretization->number_of_nodes() * 
                    energy_discretization->number_of_groups() * 
                    angular_discretization->number_of_moments(),
                    spatial_discretization->number_of_cells() * 
                    spatial_discretization->number_of_nodes() * 
                    energy_discretization->number_of_groups() * 
                    angular_discretization->number_of_moments()),
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
    Check(x.size() == column_size());
    
    switch(scattering_type_)
    {
    case FULL:
        apply_full(x);
        break;
    case COHERENT:
        apply_coherent(x);
        break;
    case INCOHERENT:
        apply_incoherent(x);
        break;
    }
    
    Check(x.size() == row_size());
}

void Scattering_Operator::
apply_incoherent(vector<double> &x)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_ordinates();

    vector<double> y(x);
    
    apply_full(y);
    
    apply_coherent(x);
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = y[k_phi] - x[k_phi];
                }
            }
        }
    }
}
