#include "RBF_Mesh.hh"

#include <vector>

#include "Check.hh"

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
