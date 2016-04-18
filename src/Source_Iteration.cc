#include "Source_Iteration.hh"

#include <cmath>

using namespace std;

Source_Iteration::
Source_Iteration(int max_iterations,
                 double tolerance,
                 shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Vector_Operator> sweeper,
                 shared_ptr<Vector_Operator> discrete_to_moment,
                 shared_ptr<Vector_Operator> moment_to_discrete,
                 shared_ptr<Vector_Operator> scattering,
                 shared_ptr<Vector_Operator> fission,
                 shared_ptr<Source_Data> source_data):
    Vector_Operator(moment_to_discrete->column_size(),
                    moment_to_discrete->column_size()),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    sweeper_(sweeper),
    discrete_to_moment_(discrete_to_moment),
    moment_to_discrete_(moment_to_discrete),
    scattering_(scattering),
    fission_(fission),
    source_data_(source_data)
{
}

/*
  Apply phi(l+1) = D(Linv(M(S phi(l)))) + D(Linv(q))
*/
void Source_Iteration::
apply(vector<double> &x)
{
    shared_ptr<Vector_Operator> Linv = sweeper_;
    shared_ptr<Vector_Operator> D = discrete_to_moment_;
    shared_ptr<Vector_Operator> M = moment_to_discrete_;
    shared_ptr<Vector_Operator> S = scattering_;
    shared_ptr<Vector_Operator> F = fission_;
    vector<double> const q_dat = source_data_->internal_source();
    
    vector<double> q(q_dat);
    
    (*D)((*Linv)(q));
    
    vector<double> x_old;
    
    for (int i = 0; i < max_iterations_; ++i)
    {
        x_old = x;

        (*S)(x);
        (*Linv)(x);
        (*D)(x);
        
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] = x[i] + q[i];
        }
        
        if (check_convergence(x, x_old))
        {
            total_iterations_ = i + 1;
            
            break;
        }
    }
}

bool Source_Iteration::
check_convergence(vector<double> const &x, 
                  vector<double> const &x_old)
{
    int number_of_cells = spatial_discretization_->number_of_cells();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    
    {
        int m = 0;
        
        for (int i = 0; i < number_of_cells; ++i)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));

                    if (abs(x[k] - x_old[k]) / (abs(x_old[k]) + tolerance_ * tolerance_) > tolerance_)
                    {
                        return false;
                    }
                }
            }
        }
    }
    
    return true;
}
