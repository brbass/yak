#ifndef DFEM_Sweep_hh
#define DFEM_Sweep_hh

#include <memory>

#include "Finite_Element_Mesh.hh"
#include "Sweep_Operator.hh"

using std::shared_ptr;

class DFEM_Sweep : public Ordinate_Sweep_Operator
{
public:
    
    DFEM_Sweep(shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Nuclear_Data> nuclear_data);
    
private:
    
    virtual void apply(vector<double> &x);
    
    unique_ptr<Finite_Element_Mesh> finite_element_mesh_;
}

#endif
