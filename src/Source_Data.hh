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
        ISOTROPIC,
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
    vector<double> const &phi_boundary() const
    {
        return phi_boundary_;
    }
    vector<double> const &psi_boundary() const
    { 
        return psi_boundary_;
    }
    
    void check_class_invariants() const;

    void update_psi_boundary(vector<double> const &psi);
    void update_phi_boundary(vector<double> const &phi);

private:

    Source_Type internal_source_type_;
    Source_Type boundary_source_type_;

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    
    vector<double> internal_source_;
    vector<double> boundary_source_;
    vector<double> alpha_;
    vector<double> phi_boundary_;
    vector<double> psi_boundary_;
};

#endif
