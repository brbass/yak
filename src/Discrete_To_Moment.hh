#ifndef Discrete_To_Moment_hh
#define Discrete_To_Moment_hh

#include "Vector_Operator.hh"

#include <vector>

using std::vector;

class Discrete_To_Moment: public Vector_Operator
{
public:
    
    Discrete_To_Moment(int number_of_cells,
                       int number_of_nodes,
                       int number_of_groups,
                       int number_of_moments,
                       int number_of_ordinates);
    
    virtual int row_size()
    {
        return number_of_cells_ * number_of_nodes_ * number_of_groups_ * number_of_moments_;
    }
    virtual int column_size()
    {
        return number_of_cells_ * number_of_nodes_ * number_of_groups_ * number_of_ordinates_;
    }
    virtual bool square()
    {
        return (row_size() == column_size());
    }
    
private:

    virtual void apply(vector<double> &x);
    
    int number_of_cells_;
    int number_of_nodes_;
    int number_of_groups_;
    int number_of_moments_;
    int number_of_ordinates_;
    
    vector<double> weights_;
    vector<double> ordinates_;
};

#endif
