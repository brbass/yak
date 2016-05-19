#ifndef GSL_Dense_Solve_hh
#define GSL_Dense_Solve_hh

#include <vector>

using std::vector;

/*
  Dense matrix solution routines from GNU Scientific Library
  
  Not optimized for matrix reuse - does not store LU or QR factorizations
*/
class GSL_Dense_Solve
{
private:
    
public:

    // Constructor
    GSL_Dense_Solve();

    // Solve linear system via LU decomposition
    void lu_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_rows);

    // Solve linear system via QR decomposition
    void qr_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_rows);

    // Solve least squares problem via QR decomposition
    void qr_lssolve(vector<double> &a_data,
                    vector<double> &b_data,
                    vector<double> &x_data,
                    unsigned number_of_rows,
                    unsigned number_of_columns);
};

#endif
