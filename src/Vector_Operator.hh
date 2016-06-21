#ifndef Vector_Operator_hh
#define Vector_Operator_hh

#include <vector>

#include "Check.hh"

using std::vector;

/*
  Pure virtual class to represent a vector operator
*/
class Vector_Operator
{
public:

    // Constructor
    Vector_Operator(int row_size,
                    int column_size);

    // Apply the operator
    vector<double> &operator()(vector<double> &x)
    {
        Check(x.size() == column_size_);
        
        apply(x);
        
        Check(x.size() == row_size_);
        
        return x;
    }

    // Output size
    virtual int row_size() const
    {
        return row_size_;
    }

    // Input size
    virtual int column_size() const
    {
        return column_size_;
    }

    // Are the input and output size the same?
    virtual bool square() const
    {
        return (row_size_ == column_size_);
    }
    
private:
    
    virtual void apply(vector<double> &x) const = 0;
    
    int row_size_;
    int column_size_;
};

#endif
