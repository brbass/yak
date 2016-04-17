#include "RBF_Mesh.hh"

#include <vector>

#include "Check.hh"

using namespace std;

RBF_Mesh::
RBF_Mesh(int dimension,
         int number_of_points,
         Geometry geometry,
         vector<double> const &positions):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_points_(number_of_points),
    point_positions_(positions)
{
    check_class_invariants();
}

void RBF_Mesh::
check_class_invariants() const
{
    Assert(point_positions_.size() == number_of_points_ * dimension_);
}
