#ifndef RBF_Solid_Mesh_hh
#define RBF_Solid_Mesh_hh

class RBF_Solid_Mesh : public virtual RBF_Mesh
{
public:

    RBF_Solid_Mesh(int dimension,
                   int number_of_points,
                   Geometry geometry,
                   Basis_Type basis_type,
                   vector<int> const &material,
                   vector<double> const &positions,
                   vector<double> const &shape_parameter);
}

#endif
