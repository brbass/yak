#ifndef Trilinos_Dense_Solve_hh
#define Trilinos_Dense_Solve_hh

#include <vector>

using std::vector;

class Trilinos_Dense_Solve
{
private:
    
public:

    Trilinos_Dense_Solve();
    
    void epetra_solve(vector<double> &a_data,
                      vector<double> &b_data,
                      vector<double> &x_data,
                      unsigned number_of_elements);
  
    void amesos_dense_solve(vector<double> &a_data,
                            vector<double> &b_data,
                            vector<double> &x_data,
                            unsigned number_of_elements);
  
    void aztec_dense_solve(vector<double> &a_data,
                           vector<double> &b_data,
                           vector<double> &x_data,
                           unsigned number_of_elements);

};

#endif
