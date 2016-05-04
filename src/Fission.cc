#include "Fission.hh"

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Fission::
Fission(shared_ptr<Spatial_Discretization> spatial_discretization,
        shared_ptr<Angular_Discretization> angular_discretization,
        shared_ptr<Energy_Discretization> energy_discretization,
        shared_ptr<Nuclear_Data> nuclear_data,
        Scattering_Type scattering_type):
    Scattering_Operator(spatial_discretization,
                        angular_discretization,
                        energy_discretization,
                        nuclear_data,
                        scattering_type)
{
}

void Fission::
apply_full(vector<double> &x)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    double angular_normalization = angular_discretization_->angular_normalization();
    vector<double> const nu = nuclear_data_->nu();
    vector<double> const sigma_f = nuclear_data_->sigma_f();
    vector<double> const chi = nuclear_data_->chi();
    
    vector<double> z(number_of_cells * number_of_nodes, 0); // fission source
    
    {
        int m = 0;
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                double sum = 0;
                
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    int k_sig = g + number_of_groups * i;
                    
                    sum += nu[k_sig] * sigma_f[k_sig] * x[k_phi];
                }
                
                int k_fis = n + number_of_nodes * i;
                
                z[k_fis] = sum;
            }
        }
    }
    
    x.assign(x.size(), 0);
    
    {
        int m = 0;
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_chi = g + number_of_groups * i;
                
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_fis = n + number_of_nodes * i;
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    
                    x[k_phi] = chi[k_chi] * z[k_fis];
                }
            }
        }
    }
}

void Fission::
apply_coherent(vector<double> &x)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    double angular_normalization = angular_discretization_->angular_normalization();
    vector<double> const nu = nuclear_data_->nu();
    vector<double> const sigma_f = nuclear_data_->sigma_f();
    vector<double> const chi = nuclear_data_->chi();
    
    vector<double> y(x);
    
    {
        int m = 0;
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_sig = g + number_of_groups * i;
                
                double cs = chi[k_sig] / angular_normalization * nu[k_sig] * sigma_f[k_sig];
                
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = cs * y[k_phi];
                }
            }
        }
    }
}

