#include "Krylov_Iteration.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziGeneralizedDavidsonSolMgr.hpp>
#include <AztecOO.h>
#include <mpi.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include "Teuchos_RCPStdSharedPtrConversions.hpp"

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Epetra_Operator_Interface.hh"
#include "Identity_Operator.hh"
#include "Nuclear_Data.hh"
#include "Preconditioner.hh"
#include "Sweep_Operator.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator_Functions.hh"
#include "XML_Functions.hh"

typedef Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> Anasazi_Eigenproblem;
typedef Anasazi::GeneralizedDavidsonSolMgr<double, Epetra_MultiVector, Epetra_Operator> Anasazi_SolverManager;
typedef Anasazi::Eigensolution<double, Epetra_MultiVector> Anasazi_Eigensolution;

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::get_shared_ptr;

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
                 shared_ptr<Vector_Operator> fission,
                 shared_ptr<Vector_Operator> preconditioner):
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
    fission_(fission),
    preconditioner_(preconditioner)
{
    if (preconditioner_ == shared_ptr<Vector_Operator>())
    {
        preconditioned_ = false;
    }
    else
    {
        preconditioned_ = true;
    }
}

void Krylov_Iteration::
solve_steady_state(vector<double> &x)
{
    shared_ptr<Vector_Operator> SI = make_shared<Source_Iterator>(*this);
    shared_ptr<Vector_Operator> FI = make_shared<Flux_Iterator>(*this);
    
    // Sweep source 
    
    vector<double> q(phi_size() + number_of_augments(), 0);
    
    print_name("Initial source iteration");
    
    if (source_data_->has_reflection())
    {
        vector<double> q_old;
        
        for (int it = 0; it < max_iterations_; ++it)
        {
            print_iteration(it);
            
            q_old = q;
            
            (*SI)(q);

            double error;
            bool converged = check_phi_convergence(q, q_old, error);
            print_error(error);
            
            if (converged)
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
        print_iteration(0);
        (*SI)(q);
        print_convergence();
        
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
                      vector<double> const &x_old,
                      double &error)
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

                    error = abs(x[k] - x_old[k]) / (abs(x_old[k]) + tolerance_ * tolerance_);
                    
                    if (error > tolerance_)
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
    shared_ptr<Epetra_Comm> comm = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    shared_ptr<Epetra_Map> map = make_shared<Epetra_Map>(phi_size() + number_of_augments(), 0, *comm);
    
    shared_ptr<Vector_Operator> FI
        = make_shared<Flux_Iterator>(*this,
                                     false); // don't include fission
    shared_ptr<Vector_Operator> NI
        = make_shared<Fission_Iterator>(*this);
    shared_ptr<Vector_Operator> EI
        = make_shared<Eigenvalue_Iterator>(*this,
                                           comm,
                                           map);
    
    // shared_ptr<Epetra_Operator> A
    //     = make_shared<Epetra_Operator_Interface>(comm,
    //                                              map,
    //                                              NI);
    // shared_ptr<Epetra_Operator> M
    //     = make_shared<Epetra_Operator_Interface>(comm,
    //                                              map,
    //                                              FI);
    shared_ptr<Epetra_Operator> Op
        = make_shared<Epetra_Operator_Interface>(comm,
                                                 map,
                                                 EI);

    int nev = 1; // number of eigenvalues requested
    int block_size = 1;
    int num_blocks = kspace_;
    
    shared_ptr<Epetra_MultiVector> init
        = make_shared<Epetra_MultiVector>(*map,
                                          block_size);
    init->PutScalar(1.0);
    
    shared_ptr<Anasazi_Eigenproblem> problem
        = make_shared<Anasazi_Eigenproblem>();
    problem->setOperator(rcp(Op));
    // problem->setA(rcp(A)); // causes problem to stall
    // problem->setM(rcp(M)); // causes problem to stall
    problem->setInitVec(rcp(init));
    problem->setNEV(nev);
    problem->setHermitian(false);
    bool worked = problem->setProblem();
    Check(worked);
    
    shared_ptr<ParameterList> params
        = get_shared_ptr(parameterList());
        // = make_shared<ParameterList>();

    // General parameters
    params->set("Maximum Iterations", max_iterations_);
    params->set("Block Size", block_size);
    params->set("Convergence Tolerance", tolerance_);
    params->set("Which", "LR");
    params->set("Maximum Restarts", 20);
    params->set("Relative Convergence Tolerance", true);
    if (solver_print_)
    {
        int verbosity = Anasazi::IterationDetails + Anasazi::TimingDetails + Anasazi::FinalSummary;
        params->set("Verbosity", verbosity);
        params->set("Output Frequency", 1);
    }

    // GeneralizedDavidson parameters
    params->set("Maximum Subspace Dimension", kspace_);
    params->set("Restart Dimension", int(ceil(kspace_ / 3)));
    params->set("Initial Guess", "User");
    
    shared_ptr<Anasazi_SolverManager> solver
        = make_shared<Anasazi_SolverManager>(rcp(problem), *params);
    
    if (solver->solve() != Anasazi::ReturnType::Converged)
    {
        cerr << "Eigenvalue solve did not converge" << endl;
    }
    
    total_iterations_ = solver->getNumIters();
    
    shared_ptr<Anasazi_Eigensolution const> solution
        = make_shared<Anasazi_Eigensolution> (problem->getSolution());
    
    k_eigenvalue = solution->Evals[0].realpart;
    
    shared_ptr<Epetra_MultiVector> sol
        = get_shared_ptr(solution->Evecs);
    
    x.resize(map->NumMyElements());
    
    sol->ExtractCopy(&x[0],
                     map->NumMyElements());
    
    x.resize(phi_size()); // remove augments
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
    
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(ki_.number_of_augments(), ki_.moment_to_discrete_);
    
    shared_ptr<Sweep_Operator> Linv =
        dynamic_pointer_cast<Sweep_Operator>(ki_.sweeper_);
    Assert(Linv);
    Linv->include_boundary_source(true);

    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    default: // Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }

    shared_ptr<Vector_Operator> E;
    if (ki_.preconditioned_)
    {
        shared_ptr<Preconditioner> Cinv
            = dynamic_pointer_cast<Preconditioner>(ki_.preconditioner_);
        Assert(Cinv);
        shared_ptr<Sweep_Operator> sweeper = Cinv->sweeper();
        sweeper->include_boundary_source(false);
        switch (sweeper->sweep_type())
        {
        case Sweep_Operator::Sweep_Type::MOMENT:
            E = Cinv;
            break;
        default: // Sweep_Operator::Sweep_Type::ORDINATE:
            E = D * Cinv * M;
            break;
        }
    }
    
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

    shared_ptr<Vector_Operator> Op;
    
    if (ki_.preconditioned_)
    {
        shared_ptr<Vector_Operator> S
            = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                              ki_.scattering_);
        shared_ptr<Vector_Operator> F
            = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                              ki_.fission_,
                                              true); // Only one of S or F needs to include the augments
        shared_ptr<Vector_Operator> I
            = make_shared<Identity_Operator>(ki_.phi_size() + ki_.number_of_augments());
        
        Op = (I + E * (S + F)) * T;
    }
    else
    {
        Op = T;
    }
    
    (*Op)(x);
}

Krylov_Iteration::Flux_Iterator::
Flux_Iterator(Krylov_Iteration const &ki,
              bool include_fission):
    Vector_Operator(ki.phi_size() + ki.number_of_augments(),
                    ki.phi_size() + ki.number_of_augments()),
    include_fission_(include_fission),
    ki_(ki)
{
    
}

void Krylov_Iteration::Flux_Iterator::
apply(vector<double> &x) const
{
    shared_ptr<Vector_Operator> D
        = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                          ki_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M
        = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                          ki_.moment_to_discrete_);

    shared_ptr<Sweep_Operator> Linv =
        dynamic_pointer_cast<Sweep_Operator>(ki_.sweeper_);
    Assert(Linv);
    Linv->include_boundary_source(false);

    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    default: // Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }

    shared_ptr<Vector_Operator> E;
    if (ki_.preconditioned_)
    {
        shared_ptr<Preconditioner> Cinv
            = dynamic_pointer_cast<Preconditioner>(ki_.preconditioner_);
        Assert(Cinv);
        shared_ptr<Sweep_Operator> sweeper = Cinv->sweeper();
        sweeper->include_boundary_source(false);
        switch (sweeper->sweep_type())
        {
        case Sweep_Operator::Sweep_Type::MOMENT:
            E = Cinv;
            break;
        default: // Sweep_Operator::Sweep_Type::ORDINATE:
            E = D * Cinv * M;
            break;
        }
    }
    
    // Add augments to operators
    
    shared_ptr<Vector_Operator> S
        = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                          ki_.scattering_);
    shared_ptr<Vector_Operator> SF;
    
    if (include_fission_)
    {
        shared_ptr<Vector_Operator> F
            = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                              ki_.fission_,
                                              true); // Only one of S or F needs to include the augments
        
        SF = S + F;
    }
    else
    {
        SF = S;
    }

    shared_ptr<Vector_Operator> I
        = make_shared<Identity_Operator>(ki_.phi_size() + ki_.number_of_augments());
    
    // Create combined operators
    
    shared_ptr<Vector_Operator> Op;

    if (ki_.preconditioned_)
    {
        Op = (I + E * SF) * (I - T * SF);
    }
    else
    {
        Op = I - T * SF;
    }
    
    (*Op)(x);
}

Krylov_Iteration::Fission_Iterator::
Fission_Iterator(Krylov_Iteration const &ki):
    Vector_Operator(ki.phi_size() + ki.number_of_augments(),
                    ki.phi_size() + ki.number_of_augments()),
    ki_(ki)
{
}

void Krylov_Iteration::Fission_Iterator::
apply(vector<double> &x) const
{
    shared_ptr<Vector_Operator> D
        = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                          ki_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M
        = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                          ki_.moment_to_discrete_);
    
    shared_ptr<Sweep_Operator> Linv =
        dynamic_pointer_cast<Sweep_Operator>(ki_.sweeper_);
    Assert(Linv);
    Linv->include_boundary_source(false);
    
    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    default: // Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }
    
    shared_ptr<Vector_Operator> F
        = make_shared<Augmented_Operator>(ki_.number_of_augments(),
                                          ki_.fission_,
                                          true); // zero out augments for fission source
    
    shared_ptr<Vector_Operator> Op
        = T * F;

    (*Op)(x);
}

Krylov_Iteration::Eigenvalue_Iterator::
Eigenvalue_Iterator(Krylov_Iteration const &ki,
                    shared_ptr<Epetra_Comm> comm,
                    shared_ptr<Epetra_Map> map):
    Vector_Operator(ki.phi_size() + ki.number_of_augments(),
                    ki.phi_size() + ki.number_of_augments()),
    ki_(ki),
    comm_(comm),
    map_(map)
{
}

void Krylov_Iteration::Eigenvalue_Iterator::
apply(vector<double> &x) const
{
    shared_ptr<Vector_Operator> FI
        = make_shared<Flux_Iterator>(ki_,
                                     false); // don't include fission
    shared_ptr<Vector_Operator> NI
        = make_shared<Fission_Iterator>(ki_);
    
    (*NI)(x);
    
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(*map_);
    shared_ptr<Epetra_Vector> rhs
        = make_shared<Epetra_Vector>(Copy,
                                     *map_,
                                     &x[0]);
    shared_ptr<Epetra_Operator> oper
        = make_shared<Epetra_Operator_Interface>(comm_,
                                                 map_,
                                                 FI);
    shared_ptr<Epetra_LinearProblem> problem
        = make_shared<Epetra_LinearProblem>(oper.get(),
                                            lhs.get(),
                                            rhs.get());
    shared_ptr<AztecOO> solver
        = make_shared<AztecOO>(*problem);
    
    solver->SetAztecOption(AZ_precond, AZ_none);
    solver->SetAztecOption(AZ_solver, AZ_gmres);
    solver->SetAztecOption(AZ_kspace, ki_.kspace_);
    solver->SetAztecOption(AZ_conv, AZ_rhs);
    if (ki_.solver_print_ && false)
    {
        solver->SetAztecOption(AZ_output, AZ_all);
    }
    else
    {
        solver->SetAztecOption(AZ_output, AZ_none);
    }
    
    lhs->PutScalar(1.0);
    
    solver->Iterate(ki_.max_iterations_, ki_.tolerance_);
    
    lhs->ExtractCopy(&x[0]);
}

