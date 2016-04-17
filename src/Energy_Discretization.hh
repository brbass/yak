#ifndef Energy_Discretization_hh
#define Energy_Discretization_hh

#include <vector>

using std::vector;

class Energy_Discretization
{
public:

    Energy_Discretization(int number_of_groups,
                          vector<double> &energy_bounds);
    Energy_Discretization(int number_of_groups);

    int number_of_groups()
    {
        return number_of_groups_;
    }
    vector<double> const &energy_bounds()
    {
        return energy_bounds_;
    }
    void check_class_invariants();
    
private:
    
    int number_of_groups_;
    vector<double> energy_bounds_;

};

#endif
