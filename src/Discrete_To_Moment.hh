#ifndef Discrete_To_Moment_hh
#define Discrete_To_Moment_hh

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::vector;

/*
  Converts the discrete angular flux to moments of the  angular flux
*/
class Discrete_To_Moment: public Vector_Operator
{
public:

    // Constructor
    Discrete_To_Moment(shared_ptr<Spatial_Discretization> spatial_discretization,
                       shared_ptr<Angular_Discretization> angular_discretization,
                       shared_ptr<Energy_Discretization> energy_discretization);
    
private:
    
    virtual void apply(vector<double> &x);

    // Output size
    int get_row_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                     shared_ptr<Angular_Discretization> angular_discretization,
                     shared_ptr<Energy_Discretization> energy_discretization);
    
    // Input size
    int get_column_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization);
    
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
