#ifndef GSL_Dense_Solve_hh
#define GSL_Dense_Solve_hh

#include <vector>

using std::vector;

class GSL_Dense_Solve
{
private:
    
public:

    GSL_Dense_Solve();
    
    void lu_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_rows);
  
    void qr_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_rows);

    void qr_lssolve(vector<double> &a_data,
                    vector<double> &b_data,
                    vector<double> &x_data,
                    unsigned number_of_rows,
                    unsigned number_of_columns);

};

#endif
