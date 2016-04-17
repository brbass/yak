#ifndef Angular_Discretization_hh
#define Angular_Discretization_hh

#include <vector>

using std::vector;

class Angular_Discretization
{
public:

    Angular_Discretization(int dimension,
                           int number_of_moments,
                           int number_of_ordinates);
    
    virtual int dimension()
    {
        return dimension_;
    }
    virtual int number_of_moments()
    {
        return number_of_moments_;
    }
    virtual int number_of_ordinates()
    {
        return number_of_ordinates_;
    }
    virtual double angular_normalization()
    {
        return angular_normalization_;
    }
    virtual vector<double> const &ordinates() const = 0;
    virtual vector<double> const &weights() const = 0;
    virtual double moment(int mom,
                          int ord);

private:
    
    int dimension_;
    int number_of_moments_;
    int number_of_ordinates_;
    
};

#endif
