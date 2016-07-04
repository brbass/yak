#include "Source_Data.hh"

#include <cmath>
#include <string>

#include "Check.hh"
#include "XML_Functions.hh"

using namespace std;

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

    number_of_augments_ = number_of_boundary_cells * number_of_nodes * number_of_ordinates * number_of_groups;
    
    calculate_partial_current();
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
    case Source_Type::FULL:
        Assert(internal_source_.size() == number_of_cells * number_of_nodes * number_of_ordinates * number_of_groups);
        break;
    case Source_Type::MOMENT:
        Assert(internal_source_.size() == number_of_cells * number_of_nodes * number_of_moments * number_of_groups);
        break;
    }
    switch(boundary_source_type_)
    {
    case Source_Type::FULL:
        Assert(boundary_source_.size() == number_of_boundary_cells * number_of_nodes * number_of_ordinates * number_of_groups);
        break;
    case Source_Type::MOMENT:
        Assert(boundary_source_.size() == number_of_boundary_cells * number_of_nodes * number_of_moments * number_of_groups);
        break;
    }
}

void Source_Data::
calculate_partial_current()
{
    int dimension = angular_discretization_->dimension();
    int number_of_boundaries = spatial_discretization_->number_of_boundary_points();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_groups = energy_discretization_->number_of_groups();
    vector<double> const weights = angular_discretization_->weights();
    vector<double> const ordinates = angular_discretization_->ordinates();
    vector<double> const boundary_normal = spatial_discretization_->boundary_normal();
    
    partial_current_.resize(number_of_boundaries * number_of_groups);

    for (int b = 0; b < number_of_boundaries; ++b)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            double current = 0;
            
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                double omega_dot_n = 0;
                
                for (int d = 0; d < dimension; ++d)
                {
                    int k_ord = d + dimension * o;
                    int k_bn = d + dimension * b;
                    
                    omega_dot_n += ordinates[k_ord] * boundary_normal[k_bn];
                }
                
                if (omega_dot_n < 0)
                {
                    int k_bs = g + number_of_groups * (o + number_of_ordinates * b);
                    current += abs(omega_dot_n) * boundary_source_[k_bs] * weights[o];
                }
            }

            int k_pc = g + number_of_groups * b;
            
            partial_current_[k_pc] = current;
        }
    }
}

bool Source_Data::
has_reflection() const
{
    for (int i = 0; i < alpha_.size(); ++i)
    {
        if (alpha_[i] != 0)
        {
            return true;
        }
    }
    
    return false;
}

vector<double> Source_Data::
total_source() const
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_boundary_cells = spatial_discretization_->number_of_boundary_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    
    vector<double> total_source(internal_source_);
    vector<int> const boundary_cells = spatial_discretization_->boundary_cells();

    switch(internal_source_type_)
    {
    case Source_Type::FULL:
        Assert(boundary_source_type_ == Source_Type::FULL);

        for (int b = 0; b < number_of_boundary_cells; ++b)
        {
            int i = boundary_cells[b];
            
            for (int n = 0; n < number_of_nodes; ++n)
            {
                int k = n + number_of_nodes * b;
                
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        int k_bs = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * b));
                        int k_ts = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                            
                        total_source[k_ts] += boundary_source_[k_bs];
                    }
                }
            }
        }
        break;
    case Source_Type::MOMENT:
        Assert(boundary_source_type_ == Source_Type::MOMENT);

        for (int b = 0; b < number_of_boundary_cells; ++b)
        {
            int i = boundary_cells[b];
            
            for (int n = 0; n < number_of_nodes; ++n)
            {
                int k = n + number_of_nodes * b;
                
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int m = 0; m < number_of_moments; ++m)
                    {
                        int k_bs = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * b));
                        int k_ts = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                            
                        total_source[k_ts] += boundary_source_[k_bs];
                    }
                }
            }
        }
        break;
    }

    return total_source;
}

void Source_Data::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("source_data");

    string internal_order;
    string boundary_order;
    switch(internal_source_type_)
    {
    case Source_Type::FULL:
        internal_order = "node-group-ordinate-cell";
        break;
    case Source_Type::MOMENT:
        internal_order = "node-group-moment-cell";
        break;
    }
    switch(boundary_source_type_)
    {
    case Source_Type::FULL:
        boundary_order = "group-ordinate-boundary_cell";
        break;
    case Source_Type::MOMENT:
        boundary_order = "group-moment-boundary_cell";
        break;
    }

    XML_Functions::append_child(source, alpha_, "alpha", "boundary_cell");
    XML_Functions::append_child(source, internal_source_, "internal_source", internal_order);
    XML_Functions::append_child(source, boundary_source_, "boundary_source", boundary_order);
}
