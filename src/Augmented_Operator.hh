#ifndef Augmented_Operator_hh
#define Augmented_Operator_hh

#include "Vector_Operator.hh"

class Augmented_Operator : public Vector_Operator
{
public:
    
    Augmented_Operator(unsigned number_of_augments,
                       shared_ptr<Vector_Operator> vector_operator);
    
private:
    
    virtual void apply(vector<double> &x);
    
    int number_of_augments_;
    shared_ptr<Vector_Operator> vector_operator_;
}

#endif
