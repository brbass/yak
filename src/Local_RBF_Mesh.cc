#include "Local_RBF_Mesh.hh"

#include "Check.hh"

Local_RBF_Mesh::
Local_RBF_Mesh(int dimension,
               int number_of_points,
               int number_of_neighbors,
               Geometry geometry,
               Basis_Type basis_type,
               vector<int> const &material,
               vector<double> const &positions,
               vector<double> const &shape_parameter):
    Spatial_Discretization(dimension,
                           geometry),
    RBF_Mesh(dimension,
             number_of_points,
             geometry,
             basis_type,
             material,
             positions,
             shape_parameter)
{
    Check(number_of_neighbors <= number_of_points);
    
    // Find nearest neighbors
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<int> local_neighbors;
        vector<double> local_distance;
        get_neighbors(i,
                      local_neighbors);
        neighbors_.push_back(local_neighbors);
    }
    
    // Initialize Trilinos

    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<int> const local_neighbors = neighbors(i);
        
        shared_ptr<Epetra_SerialDenseMatrix> mat = make_shared<Epetra_SerialDenseMatrix>(number_of_neighbors_, number_of_neighbors_);
        
        for (int j = 0; j < number_of_neighbors_; ++j)
        {
            int j1 = local_neighbors[j];
                
            shared_ptr<RBF> const equation_rbf = basis_functions(j1);
            vector<double> const equation_position = equation_rbf->position();
            
            for (int k = 0; k < number_of_neighbors_; ++k)
            {
                int k1 = local_neighbors[k];
                
                shared_ptr<RBF> const basis_rbf = basis_functions(k1);
                
                (*mat)(j, k) = basis_rbf->basis(equation_position);
            }
        }

        shared_ptr<Epetra_SerialDenseSolver> solver = make_shared<Epetra_SerialDenseSolver();
        
        solver.SetMatrix(*mat);
        solver.Factor();
        
        matrices_.push_back(mat);
        solvers_.push_back(solver);
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
    Check(b_data.size() == number_of_neighbors_);
    
    shared_ptr<Epetra_SerialDenseSolver> solver = solvers_[point];
    
    Epetra_SerialDenseVector x(number_of_neighbors_);
    Epetra_SerialDenseVector b(View, &data[0], number_of_neighbors_);
    
    solver->SetVectors(x, b);
    solver->Solve();
    
    data.assign(solver->X(), solver->X() + number_of_neighbors_);
}

/*
  Find the nearest neighbors for each point
*/
void Local_RBF_Mesh::
get_neighbors(int point,
              vector<int> &local_neighbors,
              vector<int> &local_distance)
{
    local_neighbors.resize(number_of_points_);

    vector<double> distance(number_of_points_);
    shared_ptr<RBF> const local_basis = basis_functions(point);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        shared_ptr<RBF> const other_basis = basis_functions(i);
        vector<double> const other_position = other_basis->position();

        local_neighbors[i] = i;
        distance[i] = local_basis->get_distance_squared(other_position);
    }
    
    partial_sort(local_neighbors.begin(),
                 local_neighbors.begin() + number_of_neighbors_,
                 local_neighbors.end(),
                 [&](int i, int j)
                 {
                     return distance[i] < distance[j];
                 });
    
    local_neighbors.resize(number_of_neigbors_);
}
