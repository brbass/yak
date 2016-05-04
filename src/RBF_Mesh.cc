#include "RBF_Mesh.hh"

#include <vector>

#include "Check.hh"
#include "XML_Child_Value.hh"

using namespace std;

RBF_Mesh::
RBF_Mesh(int dimension,
         int number_of_points,
         Geometry geometry,
         vector<int> const &material,
         vector<double> const &positions):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_points_(number_of_points),
    material_(material),
    point_positions_(positions)
{
    boundary_points_.push_back(0);
    boundary_points_.push_back(number_of_points_ - 1);
    boundary_nodes_.assign(2, true);
    
    for (int i = 1; i < number_of_points_ - 1; ++i)
    {
        internal_points_.push_back(i);
    }
    
    check_class_invariants();
}

void RBF_Mesh::
check_class_invariants() const
{
    Assert(boundary_points_.size() == 2);
    Assert(internal_points_.size() == number_of_points_ - 2);
    Assert(material_.size() == number_of_points_);
    Assert(point_positions_.size() == number_of_points_ * dimension_);
}

void RBF_Mesh::
output(pugi::xml_node &output_node)
{
    pugi::xml_node rbf = output_node.append_child("rbf_mesh");
    
    append_child(rbf, dimension_, "dimension");
    append_child(rbf, number_of_points_, "number_of_points");
    append_child(rbf, boundary_points_, "boundary_points", "point");
    append_child(rbf, internal_points_, "internal_points", "point");
    append_child(rbf, material_, "material", "point");
    append_child(rbf, point_positions_, "point_positions", "dimension-point");
}
