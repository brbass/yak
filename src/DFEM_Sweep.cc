#include "DFEM_Sweep.hh"

#include "Check.hh"

using namespace std;

DFEM_Sweep::
DFEM_Sweep(shared_ptr<Spatial_Discretization> spatial_discretization,
           shared_ptr<Angular_Discretization> angular_discretization,
           shared_ptr<Energy_Discretization> energy_discretization,
           shared_ptr<Nuclear_Data> nuclear_data,
           shared_ptr<Source_Data> source_data):
    Vector_Operator(spatial_discretization->number_of_cells() * 
                    spatial_discretization->number_of_nodes() * 
                    energy_discretization->number_of_groups() * 
                    angular_discretization->number_of_ordinates(),
                    spatial_discretization->number_of_cells() * 
                    spatial_discretization->number_of_nodes() * 
                    energy_discretization->number_of_groups() * 
                    angular_discretization->number_of_ordinates()),
    Ordinate_Sweep_Operator(spatial_discretization,
                            angular_discretization,
                            energy_discretization,
                            nuclear_data,
                            source_data)
{
    finite_element_mesh_ = dynamic_pointer_cast<Finite_Element_Mesh>(spatial_discretization);
}

void DFEM_Sweep::
apply(vector<double> &x)
{
    switch(spatial_discretization_->geometry())
    {
    case Spatial_Discretization::SLAB:
        sweep_slab(x);
        break;
    default:
        AssertMsg(false, "Sweep type not implemented");
        break;
    }
}

void DFEM_Sweep::
sweep_slab(vector<double> &x)
{
    vector<double> y(x);
    
    int number_of_cells = finite_element_mesh_->number_of_cells();
    int number_of_elements = finite_element_mesh_->number_of_elements();
    int number_of_nodes = finite_element_mesh_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_boundaries = finite_element_mesh_->number_of_boundary_cells();
    vector<double> const ordinates = angular_discretization_->ordinates();
    vector<double> const sigma_t = nuclear_data_->sigma_t();
    vector<double> const boundary_source = source_data_->boundary_source();
    vector<double> const alpha = source_data_->alpha();
    vector<double> const q = source_data_->internal_source();
    vector<double> const psi_boundary = source_data_->psi_boundary();
    vector<double> psi_inc(number_of_groups * number_of_ordinates);
    
    // boundary value at x = 0
    {
        int b = 0;
        int i = 0; 
        int n = 0;
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                int o1 = number_of_ordinates - o - 1;
                int k_bs = g + number_of_groups * (o + number_of_ordinates * b);
                int k_inc = g + number_of_groups * o;
                int k_ref = n + number_of_nodes * (g + number_of_groups * (o1 + number_of_ordinates * b));
                
                psi_inc[k_inc] = boundary_source[k_bs] + alpha[b] * psi_boundary[k_ref];
            }
        }
    }
    
    // sweep right
    for (int i = 0; i < number_of_elements; ++i)
    {
        Finite_Element const element = finite_element_mesh_->elements(i);
        
        double const length = element.length()[0];
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k_sig = g + number_of_groups * i;
            
            for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                int k_psi1 = 0 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                int k_psi2 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                int k_inc = g + number_of_groups * o;
                
                double a1 = ordinates[o] + sigma_t[k_sig] * length;
                double a2 = ordinates[o];
                double a3 = -ordinates[o];
                double a4 = ordinates[o] + sigma_t[k_sig] * length;
                double s1 = length * y[k_psi1] + 2 * ordinates[o] * psi_inc[k_inc];
                double s2 = length * y[k_psi2];
                
                x[k_psi1] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                x[k_psi2] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                
                psi_inc[k_inc] = x[k_psi2];
            }
        }
    }

    // boundary value at x = X
    {
        int b = 1;
        int i = number_of_cells - 1; 
        int n = 1;
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int o = 0; o < number_of_ordinates / 2; ++o)
            {
                int o1 = number_of_ordinates - o - 1;
                int k_bs = g + number_of_groups * (o + number_of_ordinates * b);
                int k_inc = g + number_of_groups * o;
                int k_ref = n + number_of_nodes * (g + number_of_groups * (o1 + number_of_ordinates * b));
                
                psi_inc[k_inc] = boundary_source[k_bs] + alpha[b] * psi_boundary[k_ref];
            }
        }
    }
    
    // sweep left
    for (int i = number_of_cells - 1; i >= 0; --i)
    {
        Finite_Element const element = finite_element_mesh_->elements(i);
        
        double const length = element.length()[0];
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k_sig = g + number_of_groups * i;
            
            for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
            {
                int k_psi1 = 0 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                int k_psi2 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                int k_inc = g + number_of_groups * o;

                double a1 = -ordinates[o] + sigma_t[k_sig] * length;
                double a2 = ordinates[o];
                double a3 = -ordinates[o];
                double a4 = -ordinates[o] + sigma_t[k_sig] * length;
                double s1 = length * y[k_psi1];
                double s2 = length * y[k_psi2] - 2 * ordinates[o] * psi_inc[k_inc];
                
                x[k_psi1] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                x[k_psi2] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                
                psi_inc[k_inc] = x[k_psi1];
            }
        }
    }
    
    source_data_->update_psi_boundary(x);
}

