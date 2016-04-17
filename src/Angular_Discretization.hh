#ifndef Angular_Discretization_hh
#define Angular_Discretization_hh

#include <vector>

using std::vector;

class Angular_Discretization
{
public:

    virtual int dimension() = 0;
    virtual int number_of_moments() = 0;
    virtual int number_of_ordinates() = 0;
    virtual vector<double> const &ordinates() = 0;
    virtual vector<double> const &weights() = 0;
    virtual double moment(int mom,
                          int ord) = 0;
};

#endif
