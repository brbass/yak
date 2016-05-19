#ifndef Ordinate_Sweep_Operator_hh
#define Ordinate_Sweep_Operator_hh

#include <memory>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

using std::shared_ptr;

/*
  Transport inverse operator
  
  \Omega \cdot \nabla \psi + \Sigma_t \psi = q
*/
class Ordinate_Sweep_Operator : virtual public Vector_Operator
{
public:

    // Constructor
    Ordinate_Sweep_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                            shared_ptr<Angular_Discretization> angular_discretization,
                            shared_ptr<Energy_Discretization> energy_discretization,
                            shared_ptr<Nuclear_Data> nuclear_data,
                            shared_ptr<Source_Data> source_data);

    // Include boundary source in sweep
    virtual void include_boundary_source(bool include_source)
    {
        include_boundary_source_ = include_source;
    }

protected:

    // Get size of input/output vector
    virtual int get_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                         shared_ptr<Angular_Discretization> angular_discretization,
                         shared_ptr<Energy_Discretization> energy_discretization,
                         shared_ptr<Source_Data> source_data);
    
    bool include_boundary_source_;
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
    shared_ptr<Source_Data> source_data_;

private:
    
    virtual void apply(vector<double> &x) = 0;
};

#endif
