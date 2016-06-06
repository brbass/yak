#include "Transport_Problem.hh"

#include "Timer.hh"
#include "XML_Child_Value.hh"

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
    Timer timer;

    timer.start();
    
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

    timer.stop();

    time_ = timer.time();
}

void Transport_Problem::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node solution = output_node.append_child("solution");

    append_child(solution, time_, "time");
    if (problem_type_ == K_EIGENVALUE)
    {
        append_child(solution, k_eigenvalue_, "k_eigenvalue");
    }
    append_child(solution, phi_, "phi", "node-group-moment-cell");
    // append_child(solution, psi_, "psi");
    
    solver_->output(output_node);
}
