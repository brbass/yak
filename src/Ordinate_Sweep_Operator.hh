#ifndef Sweep_Operator_hh
#define Sweep_Operator_hh

#include "Vector_Operator.hh"

class Sweep_Operator : public Vector_Operator
{
public:

    Sweep_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                   shared_ptr<Angular_Discretization> angular_discretization,
                   shared_ptr<Energy_Discretization> energy_discretization,
                   shared_ptr<Nuclear_Data> nuclear_data);
    
private:
    
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
}

#endif
