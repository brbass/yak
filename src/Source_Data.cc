#include "Source_Data.hh"

#include "Check.hh"

Source_Data::
Source_Data(Source_Type internal_source_type,
            Source_Type boundary_source_type,
            shared_ptr<Spatial_Discretization> spatial_discretization,
            shared_ptr<Angular_Discretization> angular_discretization,
            shared_ptr<Energy_Discretization> energy_discretization,
            vector<double> const &internal_source,
            vector<double> const &boundary_source,
            vector<double> const &alpha):
    internal_source_type_(internal_source_type),
    boundary_source_type_(boundary_source_type),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    internal_source_(internal_source),
    boundary_source_(boundary_source),
    alpha_(alpha)
    
{
    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    
    switch(boundary_source_type_)
    {
    case FULL:
        psi_boundary_.resize(number_of_boundary_cells * number_of_nodes * number_of_ordinates * number_of_groups);
        break;
    case MOMENT:
        phi_boundary_.resize(number_of_boundary_cells * number_of_nodes * number_of_moments * number_of_groups);
        break;
    }
    
    check_class_invariants();
}

void Source_Data::
check_class_invariants() const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();

    Assert(alpha_.size() == number_of_boundary_cells);
    switch(internal_source_type_)
    {
    case FULL:
        Assert(internal_source_.size() == number_of_cells * number_of_nodes * number_of_ordinates * number_of_groups);
        break;
    case MOMENT:
        Assert(internal_source_.size() == number_of_cells * number_of_nodes * number_of_moments * number_of_groups);
        break;
    }
    switch(boundary_source_type_)
    {
    case FULL:
        Assert(boundary_source_.size() == number_of_boundary_cells * number_of_ordinates * number_of_groups);
        Assert(psi_boundary_.size() == number_of_boundary_cells * number_of_nodes * number_of_ordinates * number_of_groups);
        break;
    case MOMENT:
        Assert(boundary_source_.size() == number_of_boundary_cells * number_of_moments * number_of_groups);
        Assert(phi_boundary_.size() == number_of_boundary_cells * number_of_nodes * number_of_moments * number_of_groups);
        break;
    }
}

void Source_Data::
update_psi_boundary(vector<double> const &psi)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    vector<int> boundary_cells = spatial_discretization_->boundary_cells();
    
    for (int b = 0; b < number_of_boundary_cells; ++b)
    {
        int i = boundary_cells[b];
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_psib = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * b));
                    int k_psi = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    
                    psi_boundary_[k_psib] = psi[k_psi];
                }
            }
        }
    }
}


void Source_Data::
update_phi_boundary(vector<double> const &phi)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_ordinates();
    vector<int> boundary_cells = spatial_discretization_->boundary_cells();
    
    for (int b = 0; b < number_of_boundary_cells; ++b)
    {
        int i = boundary_cells[b];
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int m = 0; m < number_of_moments; ++m)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phib = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * b));
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    phi_boundary_[k_phib] = phi[k_phi];
                }
            }
        }
    }
}

