#include "Nuclear_Data.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "XML_Functions.hh"

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

void Nuclear_Data::
check_class_invariants() const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_transition_cells = spatial_discretization_->number_of_transition_points();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    int cells_plus_transition = number_of_cells + number_of_transition_cells;
    
    Assert(sigma_t_.size() == cells_plus_transition * number_of_groups);
    Assert(sigma_s_.size() == cells_plus_transition * number_of_moments * number_of_groups * number_of_groups);
    Assert(nu_.size() == cells_plus_transition * number_of_groups);
    Assert(sigma_f_.size() == cells_plus_transition * number_of_groups);
    Assert(chi_.size() == cells_plus_transition * number_of_groups);
}

void Nuclear_Data::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node nuclear = output_node.append_child("nuclear_data");
    
    XML_Functions::append_child(nuclear, sigma_t_, "sigma_t", "group-cell");
    XML_Functions::append_child(nuclear, sigma_s_, "sigma_s", "group_from-group_to-moment-cell");
    XML_Functions::append_child(nuclear, nu_, "nu", "group-cell");
    XML_Functions::append_child(nuclear, sigma_f_, "sigma_f", "group-cell");
    XML_Functions::append_child(nuclear, chi_, "chi", "group-cell");
}
