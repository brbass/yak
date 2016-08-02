#include "Matrix_RBF_Sweep.hh"

#include <limits>

#include <Amesos.h>
#include <AztecOO.h>
#include <mpi.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include "Ifpack_Amesos.h"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Local_RBF_Diffusion.hh"
#include "Local_RBF_Mesh.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Matrix_RBF_Sweep::
Matrix_RBF_Sweep(shared_ptr<Spatial_Discretization> spatial_discretization,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                shared_ptr<Nuclear_Data> nuclear_data,
                shared_ptr<Source_Data> source_data,
                 Solver_Type solver_type):
    Sweep_Operator(Sweep_Type::ORDINATE,
                   spatial_discretization,
                   angular_discretization,
                   energy_discretization,
                   nuclear_data,
                   source_data),
    solver_type_(solver_type)
{
    reflection_tolerance_ = 1000 * numeric_limits<double>::epsilon();
    rbf_mesh_ = dynamic_pointer_cast<Local_RBF_Mesh>(spatial_discretization);
    Assert(rbf_mesh_);
    
    initialize_trilinos();
}

void Matrix_RBF_Sweep::
apply(vector<double> &x) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_boundary_points = rbf_mesh_->number_of_boundary_points();
    int number_of_internal_points = rbf_mesh_->number_of_internal_points();
    int number_of_transition_points = rbf_mesh_->number_of_transition_points();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<int> const boundary_points = rbf_mesh_->boundary_cells();
    vector<double> const boundary_normal = rbf_mesh_->boundary_normal();
    vector<int> const internal_points = rbf_mesh_->internal_cells();
    vector<int> const transition_points = rbf_mesh_->transition_cells();
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
                    set_boundary_rhs(b, i, o, g, x);
                }
                else
                {
                    set_internal_rhs(i, o, g, x);
                }
            }
            
            for (int p = 0; p < number_of_internal_points; ++p)
            {
                int i = internal_points[p];
                
                set_internal_rhs(i, o, g, x); // transition same as internal for RHS
            }
            
            for (int p = 0; p < number_of_transition_points; ++p)
            {
                int i = transition_points[p];
                
                set_internal_rhs(i, o, g, x);
            }
            
            int k = g + number_of_groups * o;
            
            // Perform matrix solve
            switch(solver_type_)
            {
            case Solver_Type::AMESOS:
                (*amesos_solver_[k])->Solve();
                
                break;
            case Solver_Type::AZTECOO:
            {
                aztec_solver_[k]->Iterate(max_iterations_, tolerance_);
                
                break;
            }
            default:
                AssertMsg(false, "Solver type not implemented");
            }
            
            // Update values for this group and ordinate
            
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

void Matrix_RBF_Sweep::
set_boundary_point(int b,
                   int i,
                   int o,
                   int g) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_augments = source_data_->number_of_augments();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    vector<int> const neighbors = rbf_mesh_->neighbors(i);
    
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
    
    int k = g + number_of_groups * o;
    
    mat_[k]->ReplaceGlobalValues(i,
                                 number_of_neighbors,
                                 &data[0],
                                 &neighbors[0]);
}

void Matrix_RBF_Sweep::
set_boundary_rhs(int b,
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
    int psi_size = row_size() - number_of_augments;
    
    vector<double> const boundary_source = source_data_->boundary_source();
    vector<double> const alpha = source_data_->alpha();
    
    vector<double> const boundary_normal = rbf_mesh_->boundary_normal();
    vector<double> local_normal(dimension);
    for (int d = 0; d < dimension; ++d)
    {
        int k_bn = d + dimension * b;
        
        local_normal[d] = boundary_normal[k_bn];
    }
    
    double rhs = 0;
    if (alpha[b] > reflection_tolerance_)
    {
        int o1 = angular_discretization_->reflect_ordinate(o,
                                                           local_normal);
        int k_ref = psi_size + g + number_of_groups * (o1 + number_of_ordinates * b);
        
        rhs = alpha[b] * x[k_ref];
    }
    
    if (include_boundary_source_)
    {
        int k_bs = g + number_of_groups * (o + number_of_ordinates * b);
        
        rhs += boundary_source[k_bs];
    }
    
    (*rhs_)[i] = rhs;
}

void Matrix_RBF_Sweep::
set_internal_point(int i,
                   int o,
                   int g) const
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
    int k_sig = g + number_of_groups * i;
    
    for (int n = 0; n < number_of_neighbors; ++n)
    {
        int j = neighbors[n];
        
        shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(j);
        
        
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
    
    int k = g + number_of_groups * o;
    
    mat_[k]->ReplaceGlobalValues(i,
                                 number_of_neighbors,
                                 &data[0],
                                 &neighbors[0]);
}

void Matrix_RBF_Sweep::
set_transition_point(int p,
                     int i,
                     int o,
                     int g) const
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<int> const neighbors = rbf_mesh_->neighbors(i);
    vector<double> const transition_normal = rbf_mesh_->transition_normal();
    vector<double> const ordinates = angular_discretization_->ordinates();
    vector<double> const sigma_t = nuclear_data_->sigma_t();
    
    shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(i);
    vector<double> const equation_position = equation_rbf->position();
    
    // Replace matrix values
    vector<double> data(number_of_neighbors, 0);
    double sigma_t_val;
    {
        double sum = 0;
        
        for (int d = 0; d < dimension; ++d)
        {
            int k_tn = d + dimension * p;
            int k_ord = d + dimension * o;
            
            sum += transition_normal[k_tn] * ordinates[k_ord];
        }
        
        if (sum > 0)
        {
            int k_sig = g + number_of_groups * i;
            
            sigma_t_val = sigma_t[k_sig];
        }
        else
        {
            int k_sig = g + number_of_groups * (number_of_points + p);
            
            sigma_t_val = sigma_t[k_sig];
        }
    }
    
    for (int n = 0; n < number_of_neighbors; ++n)
    {
        int j = neighbors[n];
        
        shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(j);
        
        double derivative = 0;
        
        for (int d = 0; d < dimension; ++d)
        {
            int k_ord = d + dimension * o;
            
            derivative += ordinates[k_ord] * basis_rbf->dbasis(d, equation_position);
        }
        
        data[n] = derivative + sigma_t_val * basis_rbf->basis(equation_position);
    }
    
    rbf_mesh_->convert_to_phi(i,
                              data);
    
    int k = g + number_of_groups * o;
    
    mat_[k]->ReplaceGlobalValues(i,
                                 number_of_neighbors,
                                 &data[0],
                                 &neighbors[0]);
}

void Matrix_RBF_Sweep::
set_internal_rhs(int i,
                 int o,
                 int g,
                 vector<double> const &x) const
{
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    
    int k_x = g + number_of_groups * (o + number_of_ordinates * i);

    (*rhs_)[i] = x[k_x];
}

void Matrix_RBF_Sweep::
initialize_trilinos()
{
    int dimension = rbf_mesh_->dimension();
    int number_of_points = rbf_mesh_->number_of_points();
    int number_of_neighbors = rbf_mesh_->number_of_neighbors();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_boundary_points = rbf_mesh_->number_of_boundary_points();
    int number_of_internal_points = rbf_mesh_->number_of_internal_points();
    int number_of_transition_points = rbf_mesh_->number_of_transition_points();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<int> const boundary_points = rbf_mesh_->boundary_cells();
    vector<double> const boundary_normal = rbf_mesh_->boundary_normal();
    vector<int> const internal_points = rbf_mesh_->internal_cells();
    vector<int> const transition_points = rbf_mesh_->transition_cells();
    vector<double> const ordinates = angular_discretization_->ordinates();
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);

    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);

    mat_.resize(number_of_groups * number_of_ordinates);
    problem_.resize(number_of_groups * number_of_ordinates);
    
    for (int g = 0; g < number_of_groups; ++g)
    {
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            int k = g + number_of_groups * o;
            
            mat_[k] = make_shared<Epetra_CrsMatrix>(Copy, *map_, number_of_neighbors, true);
            
            for (int i = 0; i < number_of_points; ++i)
            {
                vector<int> const neighbors = rbf_mesh_->neighbors(i);
                vector<double> ones(number_of_neighbors, 1);
                
                mat_[k]->InsertGlobalValues(i, number_of_neighbors, &ones[0], &neighbors[0]);
            }
            mat_[k]->FillComplete();
            mat_[k]->OptimizeStorage();
            
            problem_[k] = make_shared<Epetra_LinearProblem>(mat_[k].get(),
                                                            lhs_.get(),
                                                            rhs_.get());
        }
    }

    switch(solver_type_)
    {
    case Solver_Type::AMESOS:
    {
        amesos_solver_.resize(number_of_groups * number_of_ordinates);

        break;
    }
    case Solver_Type::AZTECOO:
    {
        aztec_solver_.resize(number_of_groups * number_of_ordinates);
        
        break;
    }
    default:
        AssertMsg(false, "Solver type not implemented");
    }
    
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
                    set_boundary_point(b, i, o, g);
                }
                else
                {
                    set_internal_point(i, o, g);
                }
            }
            for (int p = 0; p < number_of_internal_points; ++p)
            {
                int i = internal_points[p];

                set_internal_point(i, o, g);
            }
            for (int p = 0; p < number_of_transition_points; ++p)
            {
                int i = transition_points[p];
                
                set_transition_point(p, i, o, g);
            }
        }
    }
    
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;
            
            switch(solver_type_)
            {
            case Solver_Type::AMESOS:
            {
                Amesos factory;
                // Serial: Klu, Lapack, Umfpack
                // Parallel: Mumps, Superludist
                amesos_solver_[k] = make_shared<Amesos_BaseSolver*>(factory.Create("Klu", *problem_[k]));
                
                if (*amesos_solver_[k] == NULL)
                {
                    AssertMsg(false, "specified solver is not available");
                }
                
                (*amesos_solver_[k])->SymbolicFactorization();
                (*amesos_solver_[k])->NumericFactorization();
                break;
            }
            case Solver_Type::AZTECOO:
            {
                aztec_solver_[k] = make_shared<AztecOO>(*problem_[k]);
                
                // aztec_solver_[k]->SetAztecOption(AZ_precond, AZ_none);
                aztec_solver_[k]->SetAztecOption(AZ_precond, AZ_dom_decomp);
                aztec_solver_[k]->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
                // aztec_solver_[k]->SetAztecOption(AZ_graph>fill, number_of_neighbors);
                aztec_solver_[k]->SetAztecOption(AZ_solver, AZ_gmres);
                aztec_solver_[k]->SetAztecOption(AZ_kspace, 100);
                // aztec_solver_[k]->SetAztecOption(AZ_output, AZ_all);
                // aztec_solver_[k]->SetAztecOption(AZ_output, AZ_last);
                aztec_solver_[k]->SetAztecOption(AZ_output, AZ_none);
                
                double condition_number;
                aztec_solver_[k]->ConstructPreconditioner(condition_number);
                
                break;
            }
            default:
                AssertMsg(false, "Solver type not implemented");
            }
        }
    }
}
