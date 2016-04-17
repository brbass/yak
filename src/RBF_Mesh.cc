#include "RBF_Mesh.hh"

#include <vector>

using namespace std;

RBF_Mesh::
RBF_Mesh(int dimension,
         int number_of_points,
         vector<double> &positions):
    dimension_(dimension),
    number_of_points_(number_of_points),
    point_positions_(positions)
{
}
