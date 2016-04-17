#include "Scattering.hh"

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Scattering::
Scattering(shared_ptr<Spatial_Discretization> spatial_discretization,
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

virtual void Scattering::
apply_full(vector<double> &x)
{
    vector<double> y(x);
    
    x.resize(row_size());
    
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    double angular_normalization = angular_discretization_->angular_normalization();
    vector<double> const sigma_s = nuclear_data_->sigma_s();
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                    
                    for (int gf = 0; gf < number_of_groups; ++gf)
                    {
                        for (int m = 0; m < number_of_moments; ++m)
                        {
                            int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                            int k_sigma = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * i));
                            
                            double p = angular_discretization->moment(m, o);
                            
                            sum += (2 * static_cast<double>(m) + 1) / angular_normalization * sigma_s[k_sigma] * p * y[k_phi];
                        }
                    }
                    
                    int k_psi = n + number_of_nodes * (gt + number_of_groups * (o + number_of_ordinates * i));
                    
                    x[k_psi] = sum;
                }
            }
        }
    }
}

virtual void Scattering::
apply_coherent(vector<double> &x)
{
    vector<double> y(x);
    
    x.resize(row_size());
    
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    double angular_normalization = angular_discretization_->angular_normalization();
    vector<double> const sigma_s = nuclear_data_->sigma_s();
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                    
                    for (int m = 0; m < number_of_moments; ++m)
                    {
                        int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                        int k_sigma = g + number_of_groups * (g + number_of_groups * (m + number_of_moments * i));
                        
                        double p = angular_discretization->moment(m, o);
                        
                        sum += (static_cast<double>(m) + 1) / angular_normalization * sigma_s[k_sigma] * p * y[k_phi];
                    }
                    
                    int k_psi = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    
                    x[k_psi] = sum;
                }
            }
        }
    }
}

