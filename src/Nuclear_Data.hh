#ifndef Nuclear_Data_hh
#define Nuclear_Data_hh

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

class Nuclear_Data
{
public:

    Nuclear_Data(shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 vector<double> const &sigma_t,
                 vector<double> const &sigma_s,
                 vector<double> const &nu,
                 vector<double> const &sigma_f,
                 vector<double> const &chi,
                 vector<double> const &boundary_source);

    vector<double> const &sigma_t() const
    {
        return sigma_t_;
    }
    vector<double> const &sigma_s() const
    {
        return sigma_s_;
    }
    vector<double> const &nu() const
    {
        return nu_;
    }
    vector<double> const &sigma_f() const
    {
        return sigma_f_;
    }
    vector<double> const &chi() const
    {
        return chi_;
    }

    void check_class_invariants() const;

private:

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    
    vector<double> sigma_t_;
    vector<double> sigma_s_;
    vector<double> nu_;
    vector<double> sigma_f_;
    vector<double> chi_;
    
    vector<double> boundary_source_;
};

#endif
