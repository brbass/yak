#ifndef Nuclear_Data_hh
#define Nuclear_Data_hh

#include <memory>
#include <vector>

#include "pugixml.hpp"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

using std::shared_ptr;
using std::vector;

/*
  Holds cross section data for the problem
*/
class Nuclear_Data
{
public:

    // Constructor
    Nuclear_Data(shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 vector<double> const &sigma_t,
                 vector<double> const &sigma_s,
                 vector<double> const &nu,
                 vector<double> const &sigma_f,
                 vector<double> const &chi);
    
    // Total cross section
    vector<double> const &sigma_t() const
    {
        return sigma_t_;
    }

    // Scattering cross section
    vector<double> const &sigma_s() const
    {
        return sigma_s_;
    }

    // Average number of fission neutrons emitted
    vector<double> const &nu() const
    {
        return nu_;
    }

    // Fission cross section
    vector<double> const &sigma_f() const
    {
        return sigma_f_;
    }

    // Fission distribution function
    vector<double> const &chi() const
    {
        return chi_;
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
    void output(pugi::xml_node &output_node) const;

private:

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    
    vector<double> sigma_t_;
    vector<double> sigma_s_;
    vector<double> nu_;
    vector<double> sigma_f_;
    vector<double> chi_;
};

#endif
