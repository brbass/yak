#include "Energy_Discretization.hh"

#include <vector>

using namespace std;

Energy_Discretization::
Energy_Discretization(int number_of_groups):
    number_of_groups_(number_of_groups),
    energy_bounds_(number_of_groups, 0)
{
}

Energy_Discretization::
Energy_Discretization(int number_of_groups,
                      vector<double> &energy_bounds):
    number_of_groups_(number_of_groups),
    energy_bounds_(energy_bounds)
{
}
