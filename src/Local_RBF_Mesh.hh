#ifndef Local_RBF_Mesh_hh
#define Local_RBF_Mesh_hh

#include "RBF_Mesh.hh"

class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseSolver;
class Epetra_SerialDenseVector;

/*
  Mesh designed for use with radial basis functions which have sparse neighbor relationships
*/
class Local_RBF_Mesh : public RBF_Mesh
{
public:

    // Constructor
    Local_RBF_Mesh(int dimension,
                   int number_of_points,
                   int number_of_boundary_points,
                   int number_of_internal_points,
                   int number_of_neighbors,
                   Geometry geometry,
                   Basis_Type basis_type,
                   vector<int> const &material,
                   vector<int> const &boundary_points,
                   vector<int> const &internal_points,
                   vector<double> const &positions,
                   vector<double> const &shape_parameter,
                   vector<double> const &boundary_normal);
    
    // Number of neighbors for each point
    virtual int number_of_neighbors() const
    {
        return number_of_neighbors_;
    }

    // Get list of nearest neighbors to a point, sorted by distance
    virtual vector<int> const &neighbors(int i) const
    {
        return neighbors_[i];
    }

    // Convert matrix row from solution for coefficient to solution for result
    virtual void convert_to_phi(int point,
                                vector<double> &b_data);
    
protected:

    // Initilize matrices for conversion to solution solve
    virtual void get_neighbors(int point,
                               vector<int> &local_neighbors);
    
    int number_of_neighbors_;
    vector<vector<int> > neighbors_;
    vector<shared_ptr<Epetra_SerialDenseMatrix> > matrices_;
    vector<shared_ptr<Epetra_SerialDenseSolver> > solvers_;
};

#endif
