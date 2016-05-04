#ifndef Vector_Operator_hh
#define Vector_Operator_hh

#include <vector>

#include "Check.hh"

using std::vector;

class Vector_Operator
{
public:

    Vector_Operator(int row_size,
                    int column_size);

    vector<double> &operator()(vector<double> &x)
    {
        Check(x.size() == column_size_);
        
        apply(x);
        
        Check(x.size() == row_size_);
        
        return x;
    }
    
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
    
    virtual void apply(vector<double> &x) = 0;
    
    int row_size_;
    int column_size_;
};

#endif
