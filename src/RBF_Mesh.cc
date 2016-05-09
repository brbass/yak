#include "RBF_Mesh.hh"

#include <vector>

#include "Check.hh"
#include "Gaussian_RBF.hh"
#include "XML_Child_Value.hh"

using namespace std;

RBF_Mesh::
RBF_Mesh(int dimension,
         int number_of_points,
         Geometry geometry,
         Basis_Type basis_type,
         vector<int> const &material,
         vector<double> const &positions,
         vector<double> const &shape_parameter):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_points_(number_of_points),
    number_of_boundary_points_(2),
    number_of_internal_points_(number_of_points - 2),
    material_(material),
    point_positions_(positions),
    shape_parameter_(shape_parameter)
{
    boundary_points_.push_back(0);
    boundary_points_.push_back(number_of_points_ - 1);
    boundary_nodes_.assign(2, true);
    
    for (int i = 1; i < number_of_points_ - 1; ++i)
    {
        internal_points_.push_back(i);
    }
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<double> position = {positions[i]};
        vector<double> shape = {shape_parameter[i]};
        
        shared_ptr<RBF> rbf;
        
        switch(basis_type)
        {
        case GAUSSIAN:
        {
            rbf = make_shared<Gaussian_RBF>(dimension,
                                            position,
                                            shape);
            break;
        }
        default:
            AssertMsg(false, "no such type of basis function");
            break;
        }

        basis_functions_.push_back(rbf);
    }

    check_class_invariants();
}

void RBF_Mesh::
check_class_invariants() const
{
    Assert(boundary_points_.size() == 2);
    Assert(internal_points_.size() == number_of_points_ - 2);
    Assert(material_.size() == number_of_points_);
}

void RBF_Mesh::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node rbf = output_node.append_child("rbf_mesh");
    
    append_child(rbf, dimension_, "dimension");
    append_child(rbf, number_of_points_, "number_of_points");
    append_child(rbf, boundary_points_, "boundary_points", "point");
    append_child(rbf, internal_points_, "internal_points", "point");
    append_child(rbf, material_, "material", "point");
    append_child(rbf, point_positions_, "point_positions", "dimension-point");
    append_child(rbf, shape_parameter_, "shape_parameter", "dimension-point");
}
