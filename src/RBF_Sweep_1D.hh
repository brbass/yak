#ifndef RBF_Sweep_1D_hh
#define RBF_Sweep_1D_hh

#include <memory>

#include "RBF_Mesh.hh"
#include "Ordinate_Sweep_Operator.hh"

using std::shared_ptr;

/*
  Inverse operator for localized RBF sweep
  
  Creates and solves the matrix directly
*/
class RBF_Sweep_1D : public Ordinate_Sweep_Operator
{
public:

    // Constructor
    RBF_Sweep_1D(shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Nuclear_Data> nuclear_data,
                 shared_ptr<Source_Data> source_data);

protected:

    shared_ptr<RBF_Mesh> rbf_mesh_;

private:
    
    virtual void apply(vector<double> &x);
    void sweep_slab(vector<double> &x);
    
};

#endif
