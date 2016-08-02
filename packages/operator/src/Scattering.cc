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

void Scattering::
calculate_sigma_s(vector<double> &sigma_s) const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    
    sigma_s.resize(number_of_cells * number_of_groups * number_of_groups * number_of_scattering_moments);
    vector<double> const sigma_s_data = nuclear_data_->sigma_s();
    
    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_transition_cells = spatial_discretization_->number_of_transition_points();
    int number_of_internal_cells = number_of_cells - number_of_boundary_cells - number_of_transition_cells;
    
    vector<int> boundary_cells = spatial_discretization_->boundary_cells();
    vector<int> internal_cells = spatial_discretization_->internal_cells();
    vector<int> transition_cells = spatial_discretization_->transition_cells();
    
    for (int l = 0; l < number_of_scattering_moments; ++l)
    {
        for (int gt = 0; gt < number_of_groups; ++gt)
        {
            for (int gf = 0; gf < number_of_groups; ++gf)
            {
                for (int c = 0; c < number_of_boundary_cells; ++c)
                {
                    int i = boundary_cells[c];
                    int k_sigma = gf + number_of_groups * (gt + number_of_groups * (l + number_of_scattering_moments * i));
                    
                    sigma_s[k_sigma] = sigma_s_data[k_sigma];
                }

                for (int c = 0; c < number_of_internal_cells; ++c)
                {
                    int i = internal_cells[c];
                    int k_sigma = gf + number_of_groups * (gt + number_of_groups * (l + number_of_scattering_moments * i));
                    
                    sigma_s[k_sigma] = sigma_s_data[k_sigma];
                }

                for (int c = 0; c < number_of_transition_cells; ++c)
                {
                    int i = transition_cells[c];
                    int k_sigma = gf + number_of_groups * (gt + number_of_groups * (l + number_of_scattering_moments * i));
                    int k_sigma_trans = gf + number_of_groups * (gt + number_of_groups * (l + number_of_scattering_moments * (number_of_cells + c)));
                    
                    // switch (number_of_scattering_moments)
                    // {
                    // case 0:
                    sigma_s[k_sigma] = 0.5 * (sigma_s_data[k_sigma] + sigma_s_data[k_sigma_trans]);
                    //     break;
                    // default:
                    // }
                }
            }
        }
    }
}

void Scattering::
apply_full(vector<double> &x) const
{
    vector<double> y(x);
    
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();

    vector<double> sigma_s;
    
    calculate_sigma_s(sigma_s);
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                        
                    for (int gf = 0; gf < number_of_groups; ++gf)
                    {
                        int l = scattering_indices[m];
                            
                        int k_phi_from = n + number_of_nodes * (gf + number_of_groups * (m + number_of_moments * i));
                        int k_sigma = gf + number_of_groups * (gt + number_of_groups * (l + number_of_scattering_moments * i));
                        
                        sum += sigma_s[k_sigma] * y[k_phi_from];
                    }
                    
                    int k_phi_to = n + number_of_nodes * (gt + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi_to] = sum;
                }
            }
        }
    }
}

void Scattering::
apply_coherent(vector<double> &x) const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    vector<double> const sigma_s = nuclear_data_->sigma_s();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int l = scattering_indices[m];
                int k_sigma = g + number_of_groups * (g + number_of_groups * (l + number_of_scattering_moments * i));
                
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                    
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = sigma_s[k_sigma] * x[k_phi];
                }
            }
        }
    }
}

