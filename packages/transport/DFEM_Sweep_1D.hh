#ifndef DFEM_Sweep_1D_hh
#define DFEM_Sweep_1D_hh

#include <memory>

#include "Sweep_Operator.hh"

class Finite_Element_Mesh;

using std::shared_ptr;

/* 
   Inverse operator for DFEM sweep
   Uses lumped mass matrices for stability
*/
class DFEM_Sweep_1D : public Sweep_Operator
{
public:

    // Constructor
    DFEM_Sweep_1D(shared_ptr<Spatial_Discretization> spatial_discretization,
                  shared_ptr<Angular_Discretization> angular_discretization,
                  shared_ptr<Energy_Discretization> energy_discretization,
                  shared_ptr<Nuclear_Data> nuclear_data,
                  shared_ptr<Source_Data> source_data);
    
protected:

    shared_ptr<Finite_Element_Mesh> finite_element_mesh_;
    
private:

    virtual void apply(vector<double> &x) const override;
    void sweep_slab(vector<double> &x) const;
    void sweep_sphere(vector<double> &x) const;
    
};

#endif
