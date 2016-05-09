#ifndef Dense_Solve_hh
#define Dense_Solve_hh

#include "GSL_Dense_Solve.hh"
#include "Trilinos_Dense_Solve.hh"

#include <vector>

using std::vector;

class Dense_Solve
{
private:
    
    unsigned size_;
    
    GSL_Dense_Solve gsl_solver;
    Trilinos_Dense_Solve trilinos_solver;
    
public:
    
    Dense_Solve(unsigned size);
    
    void solve(vector<double> &a_data, 
               vector<double> &b_data,
               vector<double> &x_data);
    
    unsigned size()
    {
        return size_;
    }
};

#endif
