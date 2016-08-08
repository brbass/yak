#include "Diffusion_Synthetic_Acceleration.hh"

#include <AztecOO.h>
#include <mpi.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Energy_Discretization.hh"
#include "Epetra_Operator_Interface.hh"
#include "Identity_Operator.hh"
#include "Nuclear_Data.hh"
#include "Preconditioner.hh"
#include "Scattering_Operator.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"
#include "Sweep_Operator.hh"
#include "Vector_Operator_Functions.hh"

using std::make_shared;

Diffusion_Synthetic_Acceleration::
Diffusion_Synthetic_Acceleration(DSA_Type dsa_type,
                                 shared_ptr<Spatial_Discretization> spatial_discretization,
                                 shared_ptr<Angular_Discretization> angular_discretization,
                                 shared_ptr<Energy_Discretization> energy_discretization,
                                 shared_ptr<Nuclear_Data> nuclear_data,
                                 shared_ptr<Source_Data> source_data,
                                 shared_ptr<Sweep_Operator> sweeper,
                                 shared_ptr<Scattering_Operator> scattering,
                                 shared_ptr<Scattering_Operator> fission):
    Preconditioner(sweeper->row_size(),
                   sweeper->column_size(),
                   spatial_discretization,
                   angular_discretization,
                   energy_discretization,
                   nuclear_data,
                   source_data,
                   sweeper),
    dsa_type_(dsa_type),
    scattering_(scattering),
    fission_(fission)
{
    phi_size_ = (spatial_discretization->number_of_points()
                 * angular_discretization->number_of_moments()
                 * energy_discretization->number_of_groups());
    switch (dsa_type_)
    {
    case DSA_Type::MULTIGROUP:
        initialize_trilinos();
        break;
    case DSA_Type::REMOVAL:
        break;
    }
}

void Diffusion_Synthetic_Acceleration::
apply(vector<double> &x) const
{
    switch (dsa_type_)
    {
    case DSA_Type::MULTIGROUP:
        apply_multigroup(x);
        break;
    case DSA_Type::REMOVAL:
        apply_removal(x);
        break;
    }
}

void Diffusion_Synthetic_Acceleration::
apply_multigroup(vector<double> &x) const
{
    sweeper_->include_boundary_source(false);
    
    (*sweeper_)(x);
    
    for (int i = 0; i < column_size(); ++i)
    {
        (*rhs_)[i] = x[i];
    }
    
    solver_->Iterate(max_iterations_, tolerance_);
    
    lhs_->ExtractCopy(&x[0]);
}

void Diffusion_Synthetic_Acceleration::
apply_removal(vector<double> &x) const
{
    AssertMsg(false, "not yet implemented");
}

void Diffusion_Synthetic_Acceleration::
initialize_trilinos()
{
    int number_of_augments = source_data_->number_of_augments();
    
    shared_ptr<Vector_Operator> I
        = make_shared<Identity_Operator>(column_size());
    shared_ptr<Vector_Operator> S
        = make_shared<Augmented_Operator>(number_of_augments,
                                          scattering_,
                                          false);
    shared_ptr<Vector_Operator> F
        = make_shared<Augmented_Operator>(number_of_augments,
                                          fission_,
                                          true);

    shared_ptr<Vector_Operator> Cinv = sweeper_;
    
    shared_ptr<Vector_Operator> oper = I - Cinv * (S + F);
    
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(column_size(), 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    oper_ = make_shared<Epetra_Operator_Interface>(comm_,
                                                   map_,
                                                   oper);
    problem_ = make_shared<Epetra_LinearProblem>(oper_.get(),
                                                 lhs_.get(),
                                                 rhs_.get());
    solver_ = make_shared<AztecOO>(*problem_);
    
    solver_->SetAztecOption(AZ_precond, AZ_none);
    solver_->SetAztecOption(AZ_solver, AZ_gmres);
    solver_->SetAztecOption(AZ_kspace, kspace_);
    solver_->SetAztecOption(AZ_conv, AZ_rhs);
    solver_->SetAztecOption(AZ_output, AZ_none);
    
    lhs_->PutScalar(0.0);
}

