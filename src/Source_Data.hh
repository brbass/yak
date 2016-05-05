#ifndef Source_Data_hh
#define Source_Data_hh

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

class Source_Data
{
public:
    
    enum Source_Type
    {
        FULL,
        MOMENT
    };

    Source_Data(Source_Type internal_source_type,
                Source_Type boundary_source_type,
                shared_ptr<Spatial_Discretization> spatial_discretization,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                vector<double> const &internal_source,
                vector<double> const &boundary_source,
                vector<double> const &alpha);

    int number_of_augments() const
    {
        return number_of_augments_;
    }
    Source_Type internal_source_type() const
    {
        return internal_source_type_;
    }
    Source_Type boundary_source_type() const
    {
        return boundary_source_type_;
    }
    vector<double> const &internal_source() const
    {
        return internal_source_;
    }
    vector<double> const &boundary_source() const
    {
        return boundary_source_;
    }
    vector<double> const &alpha() const
    {
        return alpha_;
    }

    bool has_reflection() const;
    vector<double> total_source() const;
    
    void check_class_invariants() const;
    void output(pugi::xml_node &output_node) const;

private:

    Source_Type internal_source_type_;
    Source_Type boundary_source_type_;

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;

    int number_of_augments_;
    vector<double> internal_source_;
    vector<double> boundary_source_;
    vector<double> alpha_;
};

#endif
