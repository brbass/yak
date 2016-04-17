#ifndef Scattering_hh
#define Scattering_hh

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Scattering_Operator.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

using std::shared_ptr;
using std::vector;

class Scattering : public Scattering_Operator
{
public:
    
    Scattering(shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Nuclear_Data> nuclear_data,
               Scattering_Type scattering_type = FULL);
    
private: 
    
    virtual void apply_full(vector<double> &x);
    virtual void apply_coherent(vector<double> &x);
};

#endif
