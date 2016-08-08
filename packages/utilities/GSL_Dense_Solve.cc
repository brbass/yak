#include "GSL_Dense_Solve.hh"

#include <algorithm>
#include <iostream>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "Check.hh"

using namespace std;

GSL_Dense_Solve::
GSL_Dense_Solve()
{
}

// solves linear system Ax=b using LU decomposition: faster than Epetra for n<300
void GSL_Dense_Solve::
lu_solve(vector<double> &a_data,
         vector<double> &b_data,
         vector<double> &x_data,
         unsigned number_of_rows)
{
    unsigned number_of_columns = number_of_rows;

    Check(a_data.size() == number_of_rows * number_of_columns);
    Check(b_data.size() == number_of_rows);
    Check(x_data.size() == number_of_columns);
    
    // assign data
    gsl_matrix_view a = gsl_matrix_view_array(&a_data[0],
                                              number_of_rows,
                                              number_of_columns);
    gsl_vector_view b = gsl_vector_view_array(&b_data[0],
                                              number_of_rows);
    gsl_vector_view x = gsl_vector_view_array(&x_data[0],
                                              number_of_columns);
    
    // solve
    int s;
    gsl_permutation *p = gsl_permutation_alloc(number_of_rows);
    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, &x.vector);
    gsl_permutation_free(p);
}

// solves linear system Ax=b using QR decomposition: slower than LU
void GSL_Dense_Solve::
qr_solve(vector<double> &a_data,
         vector<double> &b_data,
         vector<double> &x_data,
         unsigned number_of_rows)
{
    unsigned number_of_columns = number_of_rows;
    
    Check(a_data.size() == number_of_rows * number_of_columns);
    Check(b_data.size() == number_of_rows);
    Check(x_data.size() == number_of_columns);
    
    unsigned t_size = number_of_rows;
    
    vector<double> t_data(t_size); // tau
    
    // assign data
    gsl_matrix_view a = gsl_matrix_view_array(&a_data[0],
                                              number_of_rows,
                                              number_of_columns);
    gsl_vector_view b = gsl_vector_view_array(&b_data[0],
                                              number_of_rows);
    gsl_vector_view x = gsl_vector_view_array(&x_data[0],
                                              number_of_columns);
    gsl_vector_view t = gsl_vector_view_array(&t_data[0],
                                              t_size);
    
    gsl_linalg_QR_decomp(&a.matrix, &t.vector);
    
    gsl_linalg_QR_solve(&a.matrix, &t.vector, &b.vector, &x.vector);
}

// solves least squares problem Ax=b using QR decomposition
void GSL_Dense_Solve::
qr_lssolve(vector<double> &a_data,
         vector<double> &b_data,
         vector<double> &x_data,
         unsigned number_of_rows,
         unsigned number_of_columns)
{
    Check(a_data.size() == number_of_rows * number_of_columns);
    Check(b_data.size() == number_of_rows);
    Check(x_data.size() == number_of_columns);
    
    unsigned t_size = min(number_of_columns, number_of_rows);

    vector<double> t_data(t_size); // tau
    vector<double> r_data(number_of_rows); // residual
    
    // assign data
    gsl_matrix_view a = gsl_matrix_view_array(&a_data[0],
                                              number_of_rows,
                                              number_of_columns);
    gsl_vector_view b = gsl_vector_view_array(&b_data[0],
                                              number_of_rows);
    gsl_vector_view x = gsl_vector_view_array(&x_data[0],
                                              number_of_columns);
    gsl_vector_view t = gsl_vector_view_array(&t_data[0],
                                              t_size);
    gsl_vector_view r = gsl_vector_view_array(&r_data[0],
                                              number_of_rows);
    
    gsl_linalg_QR_decomp(&a.matrix, &t.vector);
    
    gsl_linalg_QR_lssolve(&a.matrix, &t.vector, &b.vector, &x.vector, &r.vector);
}
