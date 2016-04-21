#include "Transport_Problem.hh"

Transport_Problem::
Transport_Problem(Problem_Type problem_type,
                  shared_ptr<Solver> solver):
    problem_type_(problem_type),
    solver_(solver)
{
}

void Transport_Problem::
solve()
{
    switch(problem_type_)
    {
    case STEADY_STATE:
        solver_->solve_steady_state(phi_);
        break;
    case K_EIGENVALUE:
        solver_->solve_k_eigenvalue(k_eigenvalue_, phi_);
        break;
    case TIME_DEPENDENT:
        solver_->solve_time_dependent(phi_);
        break;
    }
}
