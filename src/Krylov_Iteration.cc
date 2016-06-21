#include "Krylov_Iteration.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <AztecOO.h>
#include <mpi.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Epetra_Operator_Interface.hh"
#include "Nuclear_Data.hh"
#include "Ordinate_Sweep_Operator.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"
#include "XML_Functions.hh"

using namespace std;

Krylov_Iteration::
Krylov_Iteration(int max_iterations,
                 int kspace,
                 int solver_print,
                 double tolerance,
                 shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Nuclear_Data> nuclear_data,
                 shared_ptr<Source_Data> source_data,
                 shared_ptr<Vector_Operator> sweeper,
                 shared_ptr<Vector_Operator> discrete_to_moment,
                 shared_ptr<Vector_Operator> moment_to_discrete,
                 shared_ptr<Vector_Operator> scattering,
                 shared_ptr<Vector_Operator> fission):
    Solver(solver_print,
           spatial_discretization,
           angular_discretization,
           energy_discretization,
           nuclear_data,
           source_data),
    max_iterations_(max_iterations),
    kspace_(kspace),
    total_iterations_(0),
    source_iterations_(0),
    tolerance_(tolerance),
    sweeper_(sweeper),
    discrete_to_moment_(discrete_to_moment),
    moment_to_discrete_(moment_to_discrete),
    scattering_(scattering),
    fission_(fission)
{
}

void Krylov_Iteration::
solve_steady_state(vector<double> &x)
{
    shared_ptr<Vector_Operator> SI = make_shared<Source_Iterator>(*this);
    shared_ptr<Vector_Operator> FI = make_shared<Flux_Iterator>(*this);
    
    // Sweep source 

    vector<double> q(phi_size() + number_of_augments(), 0);
    
    if (source_data_->has_reflection())
    {
        vector<double> q_old;
        
        print_name("Initial source iteration");
        
        for (int it = 0; it < max_iterations_; ++it)
        {
            print_iteration(it);

            q_old = q;
            
            (*SI)(q);
            
            if (check_phi_convergence(q, q_old))
            {
                source_iterations_ = it + 1;

                print_convergence();
                
                break;
            }
        }
        for (int i = phi_size(); i < phi_size() + number_of_augments(); ++i)
        {
            q[i] = 0;
        }        
        if (source_iterations_ == max_iterations_)
        {
            print_failure();
        }

    }
    else
    {
        (*SI)(q);
        
        source_iterations_ = 1;
    }

    x.resize(phi_size() + number_of_augments(), 0);
  
    shared_ptr<Epetra_Comm> comm = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    shared_ptr<Epetra_Map> map = make_shared<Epetra_Map>(phi_size() + number_of_augments(), 0, *comm);
    shared_ptr<Epetra_Vector> lhs = make_shared<Epetra_Vector>(*map);
    shared_ptr<Epetra_Vector> rhs = make_shared<Epetra_Vector>(Copy, *map, &q[0]);
    shared_ptr<Epetra_Operator> oper = make_shared<Epetra_Operator_Interface>(comm,
                                                                              map,
                                                                              FI);
    shared_ptr<Epetra_LinearProblem> problem = make_shared<Epetra_LinearProblem>(oper.get(),
                                                                                 lhs.get(),
                                                                                 rhs.get());
    shared_ptr<AztecOO> solver = make_shared<AztecOO>(*problem);
    
    solver->SetAztecOption(AZ_precond, AZ_none);
    // solver->SetAztecOption(AZ_precond, AZ_Jacobi);
    solver->SetAztecOption(AZ_solver, AZ_gmres);
    solver->SetAztecOption(AZ_kspace, kspace_);
    solver->SetAztecOption(AZ_conv, AZ_rhs);
    if (solver_print_)
    {
        solver->SetAztecOption(AZ_output, AZ_all);
    }
    else
    {
        solver->SetAztecOption(AZ_output, AZ_none);
    }
    
    lhs->PutScalar(1.0);
    
    solver->Iterate(max_iterations_, tolerance_);
    total_iterations_ = solver->NumIters();
   
    lhs->ExtractCopy(&x[0]);
    
    x.resize(phi_size()); // remove augments
}

bool Krylov_Iteration::
check_phi_convergence(vector<double> const &x, 
                      vector<double> const &x_old)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    
    {
        int m = 0;
        
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    if (abs(x[k] - x_old[k]) / (abs(x_old[k]) + tolerance_ * tolerance_) > tolerance_)
                    {
                        return false;
                    }
                }
            }
        }
    }
    
    return true;
}

void Krylov_Iteration::
solve_k_eigenvalue(double &k_eigenvalue, 
                   vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Krylov_Iteration::
solve_time_dependent(vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Krylov_Iteration::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("krylov_iteration");
    
    XML_Functions::append_child(source, max_iterations_, "max_iterations");
    XML_Functions::append_child(source, source_iterations_, "source_iterations");
    XML_Functions::append_child(source, total_iterations_, "total_iterations");
    XML_Functions::append_child(source, tolerance_, "tolerance");
    
    nuclear_data_->output(output_node);
    source_data_->output(output_node);
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}

Krylov_Iteration::Source_Iterator::
Source_Iterator(Krylov_Iteration const &ki):
    Vector_Operator(ki.phi_size() + ki.number_of_augments(),
                    ki.phi_size() + ki.number_of_augments()),
    ki_(ki)
{
}

void Krylov_Iteration::Source_Iterator::
apply(vector<double> &x) const
{
    vector<double> const internal_source = ki_.source_data_->internal_source();
    
    shared_ptr<Ordinate_Sweep_Operator> Linv = dynamic_pointer_cast<Ordinate_Sweep_Operator>(ki_.sweeper_);
    Assert(Linv);
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.moment_to_discrete_);
    
    Linv->include_boundary_source(true);
    
    if(ki_.source_data_->internal_source_type() == Source_Data::Source_Type::FULL)
    {
        vector<double> q(internal_source);
        q.resize(q.size() + ki_.number_of_augments());
        
        (*D)(q);
        for (int i = 0; i < ki_.phi_size(); ++i)
        {
            x[i] = q[i];
        }
    }
    else
    {
        for (int i = 0; i < ki_.phi_size(); ++i)
        {
            x[i] = internal_source[i];
        }
    }
    
    (*M)(x);
    (*Linv)(x);
    (*D)(x);
}

Krylov_Iteration::Flux_Iterator::
Flux_Iterator(Krylov_Iteration const &ki):
    Vector_Operator(ki.phi_size() + ki.number_of_augments(),
                    ki.phi_size() + ki.number_of_augments()),
    ki_(ki)
{
    
}

void Krylov_Iteration::Flux_Iterator::
apply(vector<double> &x) const
{
    shared_ptr<Ordinate_Sweep_Operator> Linv = dynamic_pointer_cast<Ordinate_Sweep_Operator>(ki_.sweeper_);
    Assert(Linv);
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.moment_to_discrete_);
    shared_ptr<Vector_Operator> S = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.scattering_);
    shared_ptr<Vector_Operator> F = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.fission_);

    Linv->include_boundary_source(false);
    
    vector<double> x1(x);

    {
        vector<double> x2(x);
        
        (*S)(x); // moment scattering source
        (*F)(x2); // moment fission source
        
        for (int i = 0; i < ki_.phi_size(); ++i)
        {
            x[i] += x2[i];
        }
    }
    
    (*M)(x); // discrete fission+scattering source
    (*Linv)(x); // solution
    (*D)(x); // moment solution
    
    for (int i = 0; i < x.size(); ++i)
    {
        x[i] = x1[i] - x[i];
    }
}
