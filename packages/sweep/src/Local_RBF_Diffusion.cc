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
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_boundary_points = rbf_mesh_->number_of_boundary_points();
    int number_of_internal_points = rbf_mesh_->number_of_internal_points();
    int number_of_augments = source_data_->number_of_augments();
    int phi_size = row_size() - number_of_augments;
    vector<int> const boundary_points = rbf_mesh_->boundary_cells();
    vector<int> const internal_points = rbf_mesh_->internal_cells();
    
    for (int g = 0; g < number_of_groups; ++g)
    {
        for (int b = 0; b < number_of_boundary_points; ++b)
        {
            int i = boundary_points[b];
                
            add_boundary_point(b, i, g, x);
        }
        for (int p = 0; p < number_of_internal_points; ++p)
        {
            int i = internal_points[p];
            
            add_internal_point(i, g, x);
        }
            
        // Perform matrix solve
        switch(solver_type_)
        {
        case Solver_Type::AMESOS:
            (*amesos_solver_)->NumericFactorization();
            (*amesos_solver_)->Solve();
        
            break;
        case Solver_Type::AZTECOO:
        {
            aztec_solver_->Iterate(max_iterations_, tolerance_);
            
            break;
        }
        default:
            AssertMsg(false, "Solver type not implemented");
        }

        // Update values for this group
        
        for (int i = 0; i < number_of_points; ++i)
        {
            for (int m = 0; m < number_of_moments; ++m)
            {
                int k_x = g + number_of_groups * (m + number_of_moments * i);
                
                if (m == 0)
                {
                    x[k_x] = (*lhs_)[i];
                }
                else
                {
                    x[k_x] = 0;
                }
            }
        }
    }
    
    // Augments remain unchanged
}

void Local_RBF_Diffusion::
add_boundary_point(int b,
                   int i,
                   int g,
                   vector<double> const &x) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_augments = source_data_->number_of_augments();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    vector<int> const neighbors = rbf_mesh_->neighbors(i);
    vector<double> const sigma_t = nuclear_data_->sigma_t();

    shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(i);
    vector<double> const equation_position = equation_rbf->position();
    
    vector<double> const boundary_source = source_data_->boundary_source();
    vector<double> const partial_current = source_data_->partial_current();
    vector<double> const alpha = source_data_->alpha();

    vector<double> const boundary_normal = rbf_mesh_->boundary_normal();
    
    int k_sig = g + number_of_groups * i;
    double diffusion_coefficient = 1. / (3. * sigma_t[k_sig]);
    
    // Replace matrix values
    vector<double> data(number_of_neighbors, 0);
    
    for (int n = 0; n < number_of_neighbors; ++n)
    {
        int j = neighbors[n];
        
        shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(j);
        
        double derivative = 0;
        
        for (int d = 0; d < dimension; ++d)
        {
            int k_bn = d + dimension * b;
            
            derivative += boundary_normal[k_bn] * basis_rbf->dbasis(d, equation_position);
        }
        
        data[n] = 0.25 * (1 - alpha[b]) * basis_rbf->basis(equation_position)
            + 0.5 * (1 + alpha[b]) * diffusion_coefficient * derivative;
    }
    
    rbf_mesh_->convert_to_phi(i,
                              data);
    
    mat_->ReplaceGlobalValues(i,
                              number_of_neighbors,
                              &data[0],
                              &neighbors[0]);
    
    // Replace RHS value

    if (include_boundary_source_)
    {
        int k_pc = g + number_of_groups * b;
        (*rhs_)[i] = partial_current[k_pc];
    }
    else
    {
        (*rhs_)[i] = 0;
    }
}

void Local_RBF_Diffusion::
add_internal_point(int i,
                   int g,
                   vector<double> const &x) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<int> const neighbors = rbf_mesh_->neighbors(i);
    vector<double> const sigma_t = nuclear_data_->sigma_t();
    // vector<double> const sigma_s = nuclear_data_->diffusion_coefficient();

    int k_sig = g + number_of_groups * i;
    
    double diffusion_coefficient = 1 / (3 * sigma_t[k_sig]);
    
    shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(i);
    vector<double> const equation_position = equation_rbf->position();
    
    // Replace matrix values
    vector<double> data(number_of_neighbors, 0);
    
    for (int n = 0; n < number_of_neighbors; ++n)
    {
        int j = neighbors[n];
        
        shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(j);
        
        int k_sig = g + number_of_groups * i;
        
        double derivative = 0;
        
        for (int d = 0; d < dimension; ++d)
        {
            derivative += basis_rbf->ddbasis(d, equation_position);
        }
        
        data[n] = -diffusion_coefficient * derivative + sigma_t[k_sig] * basis_rbf->basis(equation_position);
    }
    
    rbf_mesh_->convert_to_phi(i,
                              data);
    
    mat_->ReplaceGlobalValues(i,
                              number_of_neighbors,
                              &data[0],
                              &neighbors[0]);
    
    // Replace RHS value
    
    int m = 0;
    int k_x = g + number_of_groups * (m + number_of_moments * i);
    
    (*rhs_)[i] = x[k_x];
}

void Local_RBF_Diffusion::
initialize_trilinos()
{
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();

    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    mat_ = make_shared<Epetra_CrsMatrix>(Copy, *map_, number_of_neighbors, true);

    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<int> const neighbors = rbf_mesh_->neighbors(i);
        vector<double> ones(number_of_neighbors, 1);
        
        mat_->InsertGlobalValues(i, number_of_neighbors, &ones[0], &neighbors[0]);
    }
    mat_->FillComplete();
    mat_->OptimizeStorage();
    
    problem_ = make_shared<Epetra_LinearProblem>(mat_.get(),
                                                 lhs_.get(),
                                                 rhs_.get());
    
    switch(solver_type_)
    {
    case Solver_Type::AMESOS:
        Amesos factory;
        // Serial: Klu, Lapack, Umfpack
        // Parallel: Mumps, Superludist
        amesos_solver_ = make_shared<Amesos_BaseSolver*>(factory.Create("Klu", *problem_));
        
        if (*amesos_solver_ == NULL)
        {
            AssertMsg(false, "specified solver is not available");
        }
        
        (*amesos_solver_)->SymbolicFactorization();
        
        break;
    case Solver_Type::AZTECOO:
        aztec_solver_ = make_shared<AztecOO>(*problem_);
        
        aztec_solver_->SetAztecOption(AZ_precond, AZ_Jacobi);
        aztec_solver_->SetAztecOption(AZ_poly_ord, 3);
        aztec_solver_->SetAztecOption(AZ_solver, AZ_gmres);
        aztec_solver_->SetAztecOption(AZ_kspace, 10);
        aztec_solver_->SetAztecOption(AZ_output, AZ_none);
        
        break;
    default:
        AssertMsg(false, "Solver type not implemented");
    }
}
