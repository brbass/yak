#include "Augmented_Operator.hh"

Augmented_Operator::
Augmented_Operator(unsigned number_of_augments,
                   shared_ptr<Vector_Operator> vector_operator,
                   bool zero_out_augments):
    Vector_Operator(vector_operator->row_size() + number_of_augments,
                    vector_operator->column_size() + number_of_augments),
    zero_out_augments_(zero_out_augments),
    number_of_augments_(number_of_augments),
    vector_operator_(vector_operator)
{
}

void Augmented_Operator::
apply(vector<double> &x) const
{
    int operator_row_size = vector_operator_->row_size();
    int operator_column_size = vector_operator_->column_size();
    
    vector<double> y(x.begin() + operator_column_size, x.end());
    x.resize(operator_column_size);
    (*vector_operator_)(x);
    
    if (zero_out_augments_)
    {
        x.resize(row_size(), 0);
    }
    else
    {
        x.insert(x.end(), y.begin(), y.end());
    }
}
