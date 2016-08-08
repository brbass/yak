#include "Dense_Solve.hh"

#include <iostream>
#include <vector>

#include "Check.hh"
#include "GSL_Dense_Solve.hh"
#include "Trilinos_Dense_Solve.hh"

using namespace std;

Dense_Solve::
Dense_Solve(unsigned size):
    size_(size)
{
    gsl_solver_ = make_shared<GSL_Dense_Solve>();
    trilinos_solver_ = make_shared<Trilinos_Dense_Solve>();
}

// solves linear system Ax=b
void Dense_Solve::
solve(vector<double> &a_data,
      vector<double> &b_data,
      vector<double> &x_data)
{
    Check(a_data.size() == size()*size());
    Check(b_data.size() == size());
    Check(x_data.size() == size());

    if (size() < 280) // fastest for small problems
    {
        gsl_solver_->lu_solve(a_data,
                              b_data,
                              x_data,
                              size());
    }
    else // fastest for large problems
    {
        trilinos_solver_->epetra_solve(a_data,
                                       b_data,
                                       x_data,
                                       size());
    }
}

