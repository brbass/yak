#ifndef RBF_hh
#define RBF_hh

#include <vector>

using std::vector;

class RBF
{
public:
    
    RBF(int number_of_dimensions,
        vector<double> const &position,
        vector<double> const &shape_parameter);
    
    virtual void check_class_invariants();

    virtual int number_of_dimensions() const
    {
        return number_of_dimensions_;
    }
    virtual double shape(int const dim) const
    {
        return shape_[dim];
    }
    virtual vector<double> const &position() const
    {
        return position_;
    }
    
    virtual double basis(vector<double> const &r) const = 0;
    virtual double dbasis(int dim,
                          vector<double> const &r) const = 0;
    virtual double ddbasis(int dim,
                           vector<double> const &r) const = 0;
    
protected:

    int number_of_dimensions_;
    vector<double> shape_;
    vector<double> position_;
    
    virtual double get_distance_squared(vector<double> const &r) const;
};

#endif
