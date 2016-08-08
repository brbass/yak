#ifndef Vector_Operator_Sum_hh
#define Vector_Operator_Sum_hh

#include "Vector_Operator.hh"

#include <memory>

using std::shared_ptr;

/* 
   Gives the sum of two vector operators,
   
   op1(x) + op2(x)
*/
class Vector_Operator_Sum : public Vector_Operator
{
public:

    Vector_Operator_Sum(shared_ptr<Vector_Operator> op1,
                        shared_ptr<Vector_Operator> op2);

private:

    virtual void apply(vector<double> &x) const;
    
    shared_ptr<Vector_Operator> op1_;
    shared_ptr<Vector_Operator> op2_;
};

#endif
