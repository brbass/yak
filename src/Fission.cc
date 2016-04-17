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

virtual void Fission::
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
                    
                    sum += chi[k_sig] * sigma_f[k_sig] * x[k_phi];
                }
                
                int k_fis = n + number_of_nodes * i;
                
                z[k_fis] = sum;
            }
        }
    }
    
    x.resize(row_size());
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k_fis = n + number_of_nodes * i;
            int k_chi = g + number_of_groups * i;
            
            double value = z[k_fis] * chi[k_chi] / angular_normalization;
            
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_psi = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    
                    x[k_psi] = value;
                }
            }
        }
    }
}

virtual void Fission::
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
    
    x.resize(row_size());
    
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
                    
                    double value = cs * phi[k_phi];
                    
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        int k_psi = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                        
                        x[k_psi] = value;
                    }
                }
            }
        }
    }
}

