#include "Energy_Discretization.hh"

#include <vector>

#include "Check.hh"

using namespace std;

Energy_Discretization::
Energy_Discretization(int number_of_groups):
    number_of_groups_(number_of_groups),
    energy_bounds_(number_of_groups + 1, 0)
{
    for (int g = 0; g < number_of_groups + 1; ++g)
    {
        energy_bounds_[g] = g + 1;
    }
    
    check_class_invariants();
}

Energy_Discretization::
Energy_Discretization(int number_of_groups,
                      vector<double> const &energy_bounds):
    number_of_groups_(number_of_groups),
    energy_bounds_(energy_bounds)
{
    check_class_invariants();
}

void Energy_Discretization::
check_class_invariants() const
{
    Assert(energy_bounds_.size() == number_of_groups_ + 1);
}
