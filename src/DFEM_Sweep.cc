#include "DFEM_Sweep.hh"

DFEM_Sweep::
DFEM_Sweep(shared_ptr<Spatial_Discretization> spatial_discretization,
           shared_ptr<Angular_Discretization> angular_discretization,
           shared_ptr<Energy_Discretization> energy_discretization,
           shared_ptr<Nuclear_Data> nuclear_data):
    Ordinate_Sweep_Operator(spatial_discretization,
                            angular_discretization,
                            energy_discretization,
                            nuclear_data)
{
    finite_element_mesh_ = dynamic_cast<Finite_Element_Mesh>(spatial_discretization);
}

virtual void DFEM_Sweep::
apply(vector<double> &x)
{
    vector<double> y(x);
    
    
}
