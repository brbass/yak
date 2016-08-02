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
calculate_cross_sections(vector<double> &chi,
                         vector<double> &nu,
                         vector<double> &sigma_f) const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_groups = energy_discretization_->number_of_groups();

    chi.resize(number_of_cells * number_of_groups);
    nu.resize(number_of_cells * number_of_groups);
    sigma_f.resize(number_of_cells * number_of_groups);
    
    vector<double> const chi_data = nuclear_data_->chi();
    vector<double> const nu_data = nuclear_data_->nu();
    vector<double> const sigma_f_data = nuclear_data_->sigma_f();

    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_transition_cells = spatial_discretization_->number_of_transition_points();
    int number_of_internal_cells = number_of_cells - number_of_boundary_cells - number_of_transition_cells;
    
    vector<int> boundary_cells = spatial_discretization_->boundary_cells();
    vector<int> internal_cells = spatial_discretization_->internal_cells();
    vector<int> transition_cells = spatial_discretization_->transition_cells();
    
    for (int g = 0; g < number_of_groups; ++g)
    {
        for (int c = 0; c < number_of_boundary_cells; ++c)
        {
            int i = boundary_cells[c];
            int k = g + number_of_groups * i;
                
            nu[k] = nu_data[k];
            sigma_f[k] = sigma_f_data[k];
            chi[k] = chi_data[k];
        }
            
        for (int c = 0; c < number_of_internal_cells; ++c)
        {
            int i = internal_cells[c];
            int k = g + number_of_groups * i;
                    
            nu[k] = nu_data[k];
            sigma_f[k] = sigma_f_data[k];
            chi[k] = chi_data[k];
        }
            
        for (int c = 0; c < number_of_transition_cells; ++c)
        {
            int i = transition_cells[c];
            int k = g + number_of_groups * i;
            int k_trans = g + number_of_groups * (number_of_cells + c);
                
            nu[k] = 0.5 * (nu_data[k] + nu_data[k_trans]);
            sigma_f[k] = 0.5 * (sigma_f_data[k] + sigma_f_data[k_trans]);
            chi[k] = 0.5 * (chi_data[k] + chi_data[k_trans]);
        }
    }
}

void Fission::
apply_full(vector<double> &x) const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    
    vector<double> chi;
    vector<double> nu;
    vector<double> sigma_f;
    
    calculate_cross_sections(chi,
                             nu,
                             sigma_f);
    
    // Calculate fission source
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

    // Assign flux back to the appropriate group and zeroth moment
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

    // Zero out other moments
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int m = 1; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));

                    x[k_phi] = 0;
                }
            }
        }
    }
}

void Fission::
apply_coherent(vector<double> &x) const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();

    vector<double> chi;
    vector<double> nu;
    vector<double> sigma_f;

    calculate_cross_sections(chi,
                             nu,
                             sigma_f);
    
    // Apply within-group fission to zeroth moment
    {
        int m = 0;
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_sig = g + number_of_groups * i;
                
                double cs = chi[k_sig] * nu[k_sig] * sigma_f[k_sig];
                
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = cs * x[k_phi];
                }
            }
        }
    }

    // Zero out other moments
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int m = 1; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));

                    x[k_phi] = 0;
                }
            }
        }
    }
}

