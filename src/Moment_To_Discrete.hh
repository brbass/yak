#ifndef Moment_To_Discrete_hh
#define Moment_To_Discrete_hh

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::vector;

class Moment_To_Discrete: public Vector_Operator
{
public:
    
    Moment_To_Discrete(shared_ptr<Spatial_Discretization> spatial_discretization,
                       shared_ptr<Angular_Discretization> angular_discretization,
                       shared_ptr<Energy_Discretization> energy_discretization);
    
    virtual int row_size()
    {
        return row_size_;
    }
    virtual int column_size()
    {
        return column_size_;
    }
    virtual bool square()
    {
        return (row_size_ == column_size_);
    }
    
private:

    virtual void apply(vector<double> &x);
    
    int row_size_;
    int column_size_;

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
