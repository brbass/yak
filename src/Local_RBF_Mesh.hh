#ifndef Local_RBF_Mesh_hh
#define Local_RBF_Mesh_hh

#include "RBF_Mesh.hh"

Local_RBF_Mesh : public RBF_Mesh
{
public:

    Local_RBF_Mesh(int dimension,
                   int number_of_points,
                   int number_of_neighbors,
                   Geometry geometry,
                   Basis_Type basis_type,
                   vector<int> const &material,
                   vector<double> const &positions,
                   vector<double> const &shape_parameter);

    vector<int> const &neighbors(int i) const
    {
        return neighbors_[i];
    }
    void convert_to_local(vector<double> &x);
    
protected:

    void initialize_trilinos();

    int number_of_neighbors_;
    vector<vector<int> > neighbors_;
    // vector<vector<double> > neighbor_distance_;
    vector<shared_ptr<Epetra_SerialDenseMatrix> > matrices_;
    vector<shared_ptr<Epetra_SerialDenseSolver> > solvers_;
}

#endif
