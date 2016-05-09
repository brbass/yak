#include "RBF_Sweep_1D.hh"

#include "Check.hh"

using namespace std;

RBF_Sweep_1D::
RBF_Sweep_1D(shared_ptr<Spatial_Discretization> spatial_discretization,
             shared_ptr<Angular_Discretization> angular_discretization,
             shared_ptr<Energy_Discretization> energy_discretization,
             shared_ptr<Nuclear_Data> nuclear_data,
             shared_ptr<Source_Data> source_data):
    Vector_Operator(get_size(spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             source_data),
                    get_size(spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             source_data)),
    Ordinate_Sweep_Operator(spatial_discretization,
                            angular_discretization,
                            energy_discretization,
                            nuclear_data,
                            source_data)
{
    rbf_mesh_ = dynamic_pointer_cast<RBF_Mesh>(spatial_discretization);
}

void RBF_Sweep_1D::
apply(vector<double> &x)
{
    switch(spatial_discretization_->geometry())
    {
    case Spatial_Discretization::SLAB:
        sweep_slab(x);
        break;
    default:
        AssertMsg(false, "Sweep type not implemented");
        break;
    }
}

void RBF_Sweep_1D::
sweep_slab(vector<double> &x)
{
}
