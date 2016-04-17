#ifndef Identity_Operator_hh
#define Identity_Operator_hh

#include <vector>

#include "Vector_Operator.hh"

using std::vector;

class Identity_Operator : public Vector_Operator
{
public:
 
    Identity_Operator(int size);
    
private:
    
    virtual void apply(vector<double> &x);
};

#endif
