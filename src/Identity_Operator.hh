#ifndef Identity_Operator_hh
#define Identity_Operator_hh

#include <vector>

#include "Vector_Operator.hh"

using std::vector;

class Identity_Operator : public Vector_Operator
{
public:
 
    Identity_Operator(int size);

    ~Identity_Operator(){};
   
    virtual int row_size()
    {
        return size_;
    }
    
    virtual int column_size()
    {
        return size_;
    }
    
    virtual bool square()
    {
        return true;
    }

private:
    
    virtual void apply_(vector<double> &x);
    
    int size_;
};

#endif
