#include "Krylov_Iteration.hh"

#include <cmath>

#include "Augmented_Operator.hh"
#include "Check.hh"
#include "Ordinate_Sweep_Operator.hh"
#include "XML_Child_Value.hh"

using namespace std;

Krylov_Iteration::
Krylov_Iteration(int max_iterations,
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
    total_iterations_(0),
    source_iterations_(0),
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
void Krylov_Iteration::
solve_steady_state(vector<double> &x)
{
    int number_of_augments = source_data_->number_of_augments();
    int phi_size = moment_to_discrete_->column_size();
    
    shared_ptr<Ordinate_Sweep_Operator> Linv = dynamic_pointer_cast<Ordinate_Sweep_Operator>(sweeper_);
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(number_of_augments, discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(number_of_augments, moment_to_discrete_);
    shared_ptr<Vector_Operator> S = make_shared<Augmented_Operator>(number_of_augments, scattering_);
    shared_ptr<Vector_Operator> F = make_shared<Augmented_Operator>(number_of_augments, fission_);

    vector<double> const qdat = source_data_->internal_source();
    vector<double> q(qdat);
    
    Linv->include_boundary_source(true);
    
    q.resize(q.size() + number_of_augments, 0);
    
    if(source_data_->internal_source_type() == Source_Data::FULL)
    {
        (*D)(q);
    }
    if (source_data_->has_reflection())
    {
        vector<double> q_old;
        vector<double> q_mom = q;
        
        for (int it = 0; it < max_iterations_; ++it)
        {
            q_old = q;
            
            for (int i = 0; i < phi_size; ++i)
            {
                q[i] = q_mom[i];
            }
            
            (*M)(q);
            (*Linv)(q);
            (*D)(q);
            
            if (check_phi_convergence(q, q_old))
            {
                source_iterations_ = it + 1;
                
                break;
            }
        }
    }
    else
    {
        (*Linv)(q);
        
        source_iterations_ = 1;
    }
    
    Linv->include_boundary_source(false);
    
    x.resize(phi_size + number_of_augments, 0);
    vector<double> x_old;
    vector<double> x1;
    
    for (int it = 0; it < max_iterations_; ++it)
    {
        x_old = x;
        x1 = x;
        
        (*S)(x); // moment scattering source
        (*F)(x1); // moment fission source
        
        for (int i = 0; i < phi_size; ++i)
        {
            x[i] += x1[i];
        }
        
        (*M)(x); // discrete fission+scattering source
        
        (*Linv)(x); // solution
        (*D)(x); // moment solution
        
        for (int i = 0; i < phi_size; ++i)
        {
            x[i] += q[i]; // add source
        }
        
        if (check_phi_convergence(x, x_old))
        {
            total_iterations_ = it + 1;
            
            break;
        }
    }
    
    x.resize(phi_size); // remove augments
}

/*
  Check convergence of pointwise relative error in scalar flux
*/

bool Krylov_Iteration::
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

void Krylov_Iteration::
solve_k_eigenvalue(double &k_eigenvalue, 
                   vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Krylov_Iteration::
solve_time_dependent(vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Krylov_Iteration::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("source_iteration");
    
    append_child(source, max_iterations_, "max_iterations");
    append_child(source, source_iterations_, "source_iterations");
    append_child(source, total_iterations_, "total_iterations");
    append_child(source, tolerance_, "tolerance");
    
    nuclear_data_->output(output_node);
    source_data_->output(output_node);
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}
