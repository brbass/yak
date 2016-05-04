#ifndef Ordinate_Sweep_Operator_hh
#define Ordinate_Sweep_Operator_hh

#include <memory>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

using std::shared_ptr;

class Ordinate_Sweep_Operator : virtual public Vector_Operator
{
public:

    Ordinate_Sweep_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                            shared_ptr<Angular_Discretization> angular_discretization,
                            shared_ptr<Energy_Discretization> energy_discretization,
                            shared_ptr<Nuclear_Data> nuclear_data);
    
    virtual void apply(vector<double> &x) = 0;
    virtual void include_boundary_source(include_source)
    {
        include_boundary_source_ = include_source;
    }

protected:

    virtual int get_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                         shared_ptr<Angular_Discretization> angular_discretization,
                         shared_ptr<Energy_Discretization> energy_discretization);
    
    bool include_boundary_source_;
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
};

#endif
