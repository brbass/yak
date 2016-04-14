#ifndef Moment_To_Discrete_hh
#define Moment_To_Discrete_hh

#include "Vector_Operator.hh"

#include <vector>

using std::vector;

class Moment_To_Discrete: public Vector_Operator
{
public:
    
    Moment_To_Discrete(int number_of_cells,
                       int number_of_nodes,
                       int number_of_groups,
                       int number_of_moments,
                       int number_of_ordinates);
    
    virtual int row_size()
    {
        return number_of_cells_ * number_of_nodes_ * number_of_groups_ * number_of_ordinates_;
    }
    virtual int column_size()
    {
        return number_of_cells_ * number_of_nodes_ * number_of_groups_ * number_of_moments_;
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
