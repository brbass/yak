#include "Nuclear_Data.hh"

#include "Check.hh"

Nuclear_Data::
Nuclear_Data(shared_ptr<Spatial_Discretization> spatial_discretization,
             shared_ptr<Angular_Discretization> angular_discretization,
             shared_ptr<Energy_Discretization> energy_discretization,
             vector<double> const &sigma_t,
             vector<double> const &sigma_s,
             vector<double> const &nu,
             vector<double> const &sigma_f,
             vector<double> const &chi):
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    sigma_t_(sigma_t),
    sigma_s_(sigma_s),
    nu_(nu),
    sigma_f_(sigma_f),
    chi_(chi)
{
    check_class_invariants();
}

void check_class_invariants() const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_moments = angular_discretization->number_of_moments();
    int number_of_groups = energy_discretization->number_of_groups();
    
    Assert(sigma_t_.size() == number_of_cells * number_of_moments);
    Assert(sigma_s_.size() == number_of_cells * number_of_moments * number_of_groups * number_of_groups);
    Assert(nu_.size() == number_of_cells * number_of_groups);
    Assert(sigma_f_.size() == number_of_cells * number_of_groups);
    Assert(chi_.size() == number_of_cells * number_of_groups);
}