#include "Discrete_To_Moment.hh"

#include "Check.hh"
#include "Gauss_Legendre.hh"
#include "Legendre_Polynomial.hh"

Discrete_To_Moment::
Discrete_To_Moment(shared_ptr<Spatial_Discretization> spatial_discretization,
                   shared_ptr<Angular_Discretization> angular_discretization,
                   shared_ptr<Energy_Discretization> energy_discretization):
    Vector_Operator(spatial_discretization_->number_of_cells() * 
                    spatial_discretization_->number_of_nodes() * 
                    energy_discretization_->number_of_groups() * 
                    angular_discretization_->number_of_moments(),
                    spatial_discretization_->number_of_cells() * 
                    spatial_discretization_->number_of_nodes() * 
                    energy_discretization_->number_of_groups() * 
                    angular_discretization_->number_of_ordinates()),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
}

void Discrete_To_Moment::
apply(vector<double> &x)
{
    Check(x.size() == column_size());
    
    vector<double> y(x);
    
    x.resize(row_size());
    
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    vector<double> const weights = angular_discretization_->weights();
    vector<double> const ordinates = angular_discretization_->ordinates();
    
    for (int i = 0; i < number_of_cells; ++i)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                for (int m = 0; m < number_of_moments; ++m)
                {
                    double sum = 0;
                    
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        int k = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                        
                        double p = angular_discretization_->moment(m, o);
                        
                        sum += weights[o] * p * y[k];
                    }
                    
                    int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k] = sum;
                }
            }
        }
    }
    
    Check(x.size() == row_size());
}
