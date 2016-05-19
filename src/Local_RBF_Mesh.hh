#ifndef Local_RBF_Mesh_hh
#define Local_RBF_Mesh_hh

#include "RBF_Mesh.hh"

/*
  Mesh designed for use with radial basis functions with sparse neighbor relationships
*/
Local_RBF_Mesh : public RBF_Mesh
{
public:

    // Constructor
    Local_RBF_Mesh(int dimension,
                   int number_of_points,
                   int number_of_neighbors,
                   Geometry geometry,
                   Basis_Type basis_type,
                   vector<int> const &material,
                   vector<double> const &positions,
                   vector<double> const &shape_parameter);

    // Get list of nearest neighbors to a point, sorted by distance
    vector<int> const &neighbors(int i) const
    {
        return neighbors_[i];
    }

    // Convert matrix row from solution for coefficient to solution for result
    void convert_to_phi(int point,
                        vector<double> &b_data);
    
protected:

    // Initilize matrices for conversion to solution solve
    void initialize_trilinos();

    int number_of_neighbors_;
    vector<vector<int> > neighbors_;
    // vector<vector<double> > neighbor_distance_;
    vector<shared_ptr<Epetra_SerialDenseMatrix> > matrices_;
    vector<shared_ptr<Epetra_SerialDenseSolver> > solvers_;
}

#endif
