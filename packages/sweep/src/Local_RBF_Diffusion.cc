#include "Local_RBF_Diffusion.hh"

#include <Amesos.h>
#include <AztecOO.h>
#include <mpi.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Local_RBF_Mesh.hh"
#include "Nuclear_Data.hh"
#include "Local_RBF_Mesh.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Local_RBF_Diffusion::
Local_RBF_Diffusion(shared_ptr<Spatial_Discretization> spatial_discretization,
                    shared_ptr<Angular_Discretization> angular_discretization,
                    shared_ptr<Energy_Discretization> energy_discretization,
                    shared_ptr<Nuclear_Data> nuclear_data,
                    shared_ptr<Source_Data> source_data,
                    Solver_Type solver_type):
    Sweep_Operator(Sweep_Type::MOMENT,
                   spatial_discretization,
                   angular_discretization,
                   energy_discretization,
                   nuclear_data,
                   source_data),
    solver_type_(solver_type)
{
    rbf_mesh_ = dynamic_pointer_cast<Local_RBF_Mesh>(spatial_discretization);
    Assert(rbf_mesh_);
    
    initialize_trilinos();
}

void Local_RBF_Diffusion::
apply(vector<double> &x) const
{
    
}

void Local_RBF_Diffusion::
initialize_trilinos()
{
    
}
