#include "Discrete_To_Moment.hh"

#include "Check.hh"
#include "Gauss_Legendre.hh"
#include "Legendre_Polynomial.hh"

Discrete_To_Moment::
Discrete_To_Moment(int number_of_cells,
                    int number_of_nodes,
                    int number_of_groups,
                    int number_of_moments,
                    int number_of_ordinates):
     number_of_cells_(number_of_cells),
     number_of_nodes_(number_of_nodes),
     number_of_groups_(number_of_groups),
     number_of_moments_(number_of_moments),
     number_of_ordinates_(number_of_ordinates)
{
    gauss_legendre_vec(number_of_ordinates, ordinates_, weights_);
}

void Discrete_To_Moment::
apply(vector<double> &x)
{
    Check(x.size() == column_size());
    
    vector<double> y(x);
    
    x.resize(row_size());
    
    for (int i = 0; i < number_of_cells_; ++i)
    {
        for (int g = 0; g < number_of_groups_; ++g)
        {
            for (int n = 0; n < number_of_nodes_; ++n)
            {
                for (int m = 0; m < number_of_moments_; ++m)
                {
                    double sum = 0;
                    
                    for (int o = 0; o < number_of_ordinates_; ++o)
                    {
                        int k = n + number_of_nodes_ * (g + number_of_groups_ * (o + number_of_ordinates_ * i));
                        
                        double p = legendre_polynomial(m, ordinates_[o]);
                        
                        sum += weights_[o] * p * y[k];
                    }
                    
                    int k = n + number_of_nodes_ * (g + number_of_groups_ * (m + number_of_moments_ * i));
                    
                    x[k] = sum;
                }
            }
        }
    }
    
    Check(x.size() == row_size());
}
