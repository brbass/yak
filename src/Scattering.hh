#ifndef Scattering_hh
#define Scattering_hh

#include <memory>
#include <vector>

#include "Scattering_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Nuclear_Data;
class Spatial_Discretization;

using std::shared_ptr;
using std::vector;

/*
  Applies scattering to a moment representation of the flux
*/
class Scattering : public Scattering_Operator
{
public:

    // Constructor
    Scattering(shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Nuclear_Data> nuclear_data,
               Scattering_Type scattering_type = Scattering_Type::FULL);
    
private: 

    // Apply within-group and out-of-group scattering
    virtual void apply_full(vector<double> &x);

    // Apply only within-group scattering
    virtual void apply_coherent(vector<double> &x);
};

#endif
