#include "Local_RBF_Sweep.hh"

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
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Local_RBF_Sweep::
Local_RBF_Sweep(shared_ptr<Spatial_Discretization> spatial_discretization,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                shared_ptr<Nuclear_Data> nuclear_data,
                shared_ptr<Source_Data> source_data,
                Solver_Type solver_type):
    Ordinate_Sweep_Operator(spatial_discretization,
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

void Local_RBF_Sweep::
apply(vector<double> &x) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_boundary_points = rbf_mesh_->number_of_boundary_points();
    int number_of_internal_points = rbf_mesh_->number_of_internal_points();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<int> const boundary_points = rbf_mesh_->boundary_cells();
    vector<double> const boundary_normal = rbf_mesh_->boundary_normal();
    vector<int> const internal_points = rbf_mesh_->internal_cells();
    vector<double> const ordinates = angular_discretization_->ordinates();
    
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int b = 0; b < number_of_boundary_points; ++b)
            {
                int i = boundary_points[b];
                
                double sum = 0;
                
                for (int d = 0; d < dimension; ++d)
                {
                    int k_bn = d + dimension * b;
                    int k_ord = d + dimension * o;
                    
                    sum += boundary_normal[k_bn] * ordinates[k_ord];
                }
                
                if (sum < 0)
                {
                    add_boundary_point(b, i, o, g, x);
                }
                else
                {
                    add_internal_point(i, o, g, x);
                }
            }
            for (int p = 0; p < number_of_internal_points; ++p)
            {
                int i = internal_points[p];

                add_internal_point(i, o, g, x);
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

            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                
                x[k_x] = (*lhs_)[i];
            }
        }
    }

    // Update augments

    for (int b = 0; b < number_of_boundary_points; ++b)
    {
        int i = boundary_points[b];
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_b = psi_size + g + number_of_groups * (o + number_of_ordinates * b);
                int k_psi = g + number_of_groups * (o + number_of_ordinates * i);
                    
                x[k_b] = x[k_psi];
            }
        }
    }
}

void Local_RBF_Sweep::
add_boundary_point(int b,
                   int i,
                   int o,
                   int g,
                   vector<double> const &x) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_augments = source_data_->number_of_augments();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    vector<int> const neighbors = rbf_mesh_->neighbors(i);

    vector<double> const boundary_source = source_data_->boundary_source();
    vector<double> const alpha = source_data_->alpha();

    vector<double> const boundary_normal = rbf_mesh_->boundary_normal();
    vector<double const>::iterator it = boundary_normal.begin() + dimension * b;
    vector<double> const local_normal(it, it + dimension);
    
    // Replace matrix values
    vector<double> data(number_of_neighbors, 0);

    bool point_found = false;
    for (int n = 0; n < number_of_neighbors; ++n)
    {
        int j = neighbors[n];
        
        if (j == i)
        {
            data[n] = 1;
            
            point_found = true;
            break;
        }
    }
    Assert(point_found);
    
    mat_->ReplaceGlobalValues(i,
                              number_of_neighbors,
                              &data[0],
                              &neighbors[0]);

    // Replace RHS value
    int psi_size = row_size() - number_of_augments;
    int o1 = angular_discretization_->reflect_ordinate(o,
                                                       local_normal);
    int k_ref = psi_size + g + number_of_groups * (o1 + number_of_ordinates * b);
    double rhs = alpha[b] * x[k_ref];

    if (include_boundary_source_)
    {
        int k_bs = g + number_of_groups * (o + number_of_ordinates * b);

        rhs += boundary_source[k_bs];
    }
    
    (*rhs_)[i] = rhs;
}

void Local_RBF_Sweep::
add_internal_point(int i,
                   int o,
                   int g,
                   vector<double> const &x) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<int> const neighbors = rbf_mesh_->neighbors(i);
    vector<double> const ordinates = angular_discretization_->ordinates();
    vector<double> const sigma_t = nuclear_data_->sigma_t();

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
            int k_ord = d + dimension * o;
            
            derivative += ordinates[k_ord] * basis_rbf->dbasis(d, equation_position);
        }
        
        data[n] = derivative + sigma_t[k_sig] * basis_rbf->basis(equation_position);
    }
    
    rbf_mesh_->convert_to_phi(i,
                              data);
    
    mat_->ReplaceGlobalValues(i,
                              number_of_neighbors,
                              &data[0],
                              &neighbors[0]);

    // Replace RHS value
    
    int k_x = g + number_of_groups * (o + number_of_ordinates * i);
    
    (*rhs_)[i] = x[k_x];
}

void Local_RBF_Sweep::
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
