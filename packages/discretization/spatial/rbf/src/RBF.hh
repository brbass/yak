#ifndef RBF_hh
#define RBF_hh

#include <vector>

using std::vector;

/*
  Pure virtual class to represent a radial basis function
*/
class RBF
{
public:

    // Constructor
    RBF(int number_of_dimensions,
        vector<double> const &position,
        vector<double> const &shape_parameter);

    // Check class invariants
    virtual void check_class_invariants();

    // Number of dimensions
    virtual int number_of_dimensions() const
    {
        return number_of_dimensions_;
    }

    // Shape parameter
    virtual double shape(int const dim) const
    {
        return shape_[dim];
    }

    // Position of RBF node
    virtual vector<double> const &position() const
    {
        return position_;
    }

    // Value of basis function at the point r
    virtual double basis(vector<double> const &r) const = 0;

    // Derivative of basis function at the point r
    virtual double dbasis(int dim,
                          vector<double> const &r) const = 0;

    // Second derivative of the basis function at the point r
    virtual double ddbasis(int dim,
                           vector<double> const &r) const = 0;

    // Return (shape \cdot (r - position))^2
    virtual double get_distance_squared(vector<double> const &r) const;

    virtual void multiply_shape_parameter(double t);
    
protected:

    int number_of_dimensions_;
    vector<double> shape_;
    vector<double> position_;
    
};

#endif
