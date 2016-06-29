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

    enum class Coefficient_Type
    {
        PHI,
        ALPHA
    };
    
    // Constructor
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
                   vector<double> const &boundary_normal,
                   vector<int> const &surface = vector<int>(),
                   vector<int> const &region = vector<int>(),
                   shared_ptr<Solid_Geometry> const solid_geometry = shared_ptr<Solid_Geometry>());
    
    // Convert matrix row from solution for coefficient to solution for result
    virtual void convert_to_phi(int point,
                                vector<double> &b_data);

    virtual Coefficient_Type coefficient_type()
    {
        return coefficient_type_;
    }
    
protected:

    Coefficient_Type coefficient_type_;
    vector<shared_ptr<Epetra_SerialDenseMatrix> > matrices_;
    vector<shared_ptr<Epetra_SerialDenseSolver> > solvers_;
};

#endif
