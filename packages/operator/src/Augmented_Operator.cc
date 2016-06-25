#include "Augmented_Operator.hh"

Augmented_Operator::
Augmented_Operator(unsigned number_of_augments,
                   shared_ptr<Vector_Operator> vector_operator):
    Vector_Operator(vector_operator->row_size() + number_of_augments,
                    vector_operator->column_size() + number_of_augments),
    number_of_augments_(number_of_augments),
    vector_operator_(vector_operator)
{
}

void Augmented_Operator::
apply(vector<double> &x) const
{
    int row_size = vector_operator_->row_size();
    int column_size = vector_operator_->column_size();
    
    vector<double> y(x.begin() + column_size, x.end());
    x.resize(column_size);
    (*vector_operator_)(x);
    x.insert(x.end(), y.begin(), y.end());
}
