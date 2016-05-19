#ifndef Dense_Solve_hh
#define Dense_Solve_hh

#include "GSL_Dense_Solve.hh"
#include "Trilinos_Dense_Solve.hh"

#include <vector>

using std::vector;

/*
  Generalized class for "dumb" dense matrix solution
  
  Automatically uses GSL or Trilinos solvers, depending on size of problem
*/
class Dense_Solve
{
public:

    // Creator
    Dense_Solve(unsigned size);

    // Solve problem Ax=b
    void solve(vector<double> &a_data, 
               vector<double> &b_data,
               vector<double> &x_data);

    // Size of square matrix
    unsigned size()
    {
        return size_;
    }

private:
    
    unsigned size_;
    
    GSL_Dense_Solve gsl_solver;
    Trilinos_Dense_Solve trilinos_solver;
};

#endif
