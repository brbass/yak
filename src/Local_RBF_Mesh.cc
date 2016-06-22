#include "Local_RBF_Mesh.hh"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_SerialDenseVector.h>

#include "Check.hh"

using namespace std;

Local_RBF_Mesh::
Local_RBF_Mesh(int dimension,
               int number_of_points,
               int number_of_boundary_points,
               int number_of_internal_points,
               int number_of_neighbors,
               double shape_multiplier,
               Geometry geometry,
               Basis_Type basis_type,
               Coefficient_Type coefficient_type,
               vector<int> const &material,
               vector<int> const &boundary_points,
               vector<int> const &internal_points,
               vector<double> const &positions,
               vector<double> const &boundary_normal):
    RBF_Mesh(dimension,
             number_of_points,
             number_of_boundary_points,
             number_of_internal_points,
             number_of_neighbors,
             shape_multiplier,
             geometry,
             basis_type,
             material,
             boundary_points,
             internal_points,
             positions,
             boundary_normal),
    coefficient_type_(coefficient_type)
{
    // Initialize Trilinos if needed

    if (coefficient_type_ == Coefficient_Type::PHI)
    {
        for (int i = 0; i < number_of_points_; ++i)
        {
            vector<int> const local_neighbors = neighbors(i);
            
            shared_ptr<Epetra_SerialDenseMatrix> mat = make_shared<Epetra_SerialDenseMatrix>(number_of_neighbors_, number_of_neighbors_);
            
            for (int j = 0; j < number_of_neighbors_; ++j)
            {
                int j1 = local_neighbors[j];
                
                shared_ptr<RBF> const equation_rbf = basis_function(j1);
                vector<double> const equation_position = equation_rbf->position();
                
                for (int k = 0; k < number_of_neighbors_; ++k)
                {
                    int k1 = local_neighbors[k];
                    
                    shared_ptr<RBF> const basis_rbf = basis_function(k1);
                    
                    (*mat)(j, k) = basis_rbf->basis(equation_position);
                }
            }
            
            shared_ptr<Epetra_SerialDenseSolver> solver = make_shared<Epetra_SerialDenseSolver>();
            
            solver->SetMatrix(*mat);
            solver->Factor();
            
            matrices_.push_back(mat);
            solvers_.push_back(solver);
        }
    }
}

/*
  Convert a row of the full problem matrix from solving for coefficients to solving 
  for the result
*/
void Local_RBF_Mesh::
convert_to_phi(int point,
               vector<double> &data)
{
    Check(data.size() == number_of_neighbors_);
    
    shared_ptr<Epetra_SerialDenseSolver> solver = solvers_[point];
    
    Epetra_SerialDenseVector x(number_of_neighbors_);
    Epetra_SerialDenseVector b(View, &data[0], number_of_neighbors_);
    
    solver->SetVectors(x, b);
    solver->Solve();
    
    // data.assign(solver->X(), solver->X() + number_of_neighbors_);
}
