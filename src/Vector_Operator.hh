#ifndef Vector_Operator_hh
#define Vector_Operator_hh

#include <vector>

using std::vector;

class Vector_Operator
{
public:

    vector<double> &operator()(vector<double> &x)
    {
        apply(x);
        
        return x;
    }
    
    virtual int row_size() = 0;
    
    virtual int column_size() = 0;
    
    virtual bool square() = 0;
    
private:
    
    virtual void apply(vector<double> &x) = 0;
};

#endif
