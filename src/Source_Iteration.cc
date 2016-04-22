#include "Source_Iteration.hh"

#include <cmath>

#include "Check.hh"
#include "XML_Child_Value.hh"

using namespace std;

Source_Iteration::
Source_Iteration(int max_iterations,
                 double tolerance,
                 shared_ptr<Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Nuclear_Data> nuclear_data,
                 shared_ptr<Source_Data> source_data,
                 shared_ptr<Vector_Operator> sweeper,
                 shared_ptr<Vector_Operator> discrete_to_moment,
                 shared_ptr<Vector_Operator> moment_to_discrete,
                 shared_ptr<Vector_Operator> scattering,
                 shared_ptr<Vector_Operator> fission):
    Solver(spatial_discretization,
           angular_discretization,
           energy_discretization,
           nuclear_data,
           source_data),
    max_iterations_(max_iterations),
    tolerance_(tolerance),
    sweeper_(sweeper),
    discrete_to_moment_(discrete_to_moment),
    moment_to_discrete_(moment_to_discrete),
    scattering_(scattering),
    fission_(fission)
{
}

/*
  Apply phi(l+1) = D(Linv(M(S phi(l)))) + D(Linv(q))
*/
void Source_Iteration::
solve_steady_state(vector<double> &x)
{
    shared_ptr<Vector_Operator> Linv = sweeper_;
    shared_ptr<Vector_Operator> D = discrete_to_moment_;
    shared_ptr<Vector_Operator> M = moment_to_discrete_;
    shared_ptr<Vector_Operator> S = scattering_;
    shared_ptr<Vector_Operator> F = fission_;
    vector<double> const q_dat = source_data_->internal_source();
    
    vector<double> q(q_dat);
    
    // if(source_data_->internal_source_type() == Source_Data::MOMENT)
    // {
    //     (*M)(q);
    // }
    
    // (*D)((*Linv)(q));
    // (*D)(q);

    x.resize(moment_to_discrete_->column_size(), 0);
    vector<double> x_old;
    vector<double> x1;

    for (int it = 0; it < max_iterations_; ++it)
    {
        x_old = x;
        x1 = x;
        
        (*S)(x); // moment scattering source
        (*F)(x1); // moment fission source
        
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] += x1[i];
        }
        
        (*M)(x); // discrete fission+scattering source
        
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] += q[i]; // discrete total source
        }
        
        (*Linv)(x); // solution
        (*D)(x); // moment solution
        
        if (check_phi_convergence(x, x_old))
        {
            total_iterations_ = it + 1;
            
            break;
        }
    }
}

/*
  Check convergence of pointwise relative error in scalar flux
*/

bool Source_Iteration::
check_phi_convergence(vector<double> const &x, 
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

void Source_Iteration::
solve_k_eigenvalue(double &k_eigenvalue, 
                   vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Source_Iteration::
solve_time_dependent(vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Source_Iteration::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("source_iteration");
    
    append_child(source, max_iterations_, "max_iterations");
    append_child(source, total_iterations_, "total_iterations");
    append_child(source, tolerance_, "tolerance");
    
    nuclear_data_->output(output_node);
    source_data_->output(output_node);
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}
