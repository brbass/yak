#include "Source_Iteration.hh"

#include <cmath>

#include "Augmented_Operator.hh"
#include "Check.hh"
#include "Ordinate_Sweep_Operator.hh"
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
void Source_Iteration::
solve_steady_state(vector<double> &x)
{
    shared_ptr<Vector_Operator> SI = make_shared<Source_Iterator>(*this);
    shared_ptr<Vector_Operator> FI = make_shared<Flux_Iterator>(*this);
    
    vector<double> q(phi_size() + number_of_augments(), 0);
    
    if (source_data_->has_reflection())
    {
        vector<double> q_old;
        
        for (int it = 0; it < max_iterations_; ++it)
        {
            q_old = q;
            
            (*SI)(q);
            
            if (check_phi_convergence(q, q_old))
            {
                source_iterations_ = it + 1;
                
                break;
            }
        }
        for (int i = phi_size(); i < phi_size() + number_of_augments(); ++i)
        {
            q[i] = 0;
        }
    }
    else
    {
        (*SI)(q);
        
        source_iterations_ = 1;
    }
    
    x.resize(phi_size() + number_of_augments(), 0);
    vector<double> x_old;
    
    for (int it = 0; it < max_iterations_; ++it)
    {
        x_old = x;
        
        (*FI)(x);
        
        for (int i = 0; i < phi_size(); ++i)
        {
            x[i] += q[i]; // add source
        }
        
        if (check_phi_convergence(x, x_old))
        {
            total_iterations_ = it + 1;
            
            break;
        }
    }
    
    x.resize(phi_size()); // remove augments
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
    append_child(source, source_iterations_, "source_iterations");
    append_child(source, total_iterations_, "total_iterations");
    append_child(source, tolerance_, "tolerance");
    
    nuclear_data_->output(output_node);
    source_data_->output(output_node);
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}

Source_Iteration::Source_Iterator::
Source_Iterator(Source_Iteration const &si):
    Vector_Operator(si.phi_size() + si.number_of_augments(),
                    si.phi_size() + si.number_of_augments()),
    si_(si)
{
}

void Source_Iteration::Source_Iterator::
apply(vector<double> &x)
{
    vector<double> const internal_source = si_.source_data_->internal_source();
    
    shared_ptr<Ordinate_Sweep_Operator> Linv = dynamic_pointer_cast<Ordinate_Sweep_Operator>(si_.sweeper_);
    Assert(Linv);
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(si_.number_of_augments(), si_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(si_.number_of_augments(), si_.moment_to_discrete_);
    
    Linv->include_boundary_source(true);
    
    if(si_.source_data_->internal_source_type() == Source_Data::FULL)
    {
        vector<double> q(internal_source);
        q.resize(q.size() + si_.number_of_augments());
        
        (*D)(q);
        for (int i = 0; i < si_.phi_size(); ++i)
        {
            x[i] = q[i];
        }
    }
    else
    {
        for (int i = 0; i < si_.phi_size(); ++i)
        {
            x[i] = internal_source[i];
        }
    }
    
    (*M)(x);
    (*Linv)(x);
    (*D)(x);
}

Source_Iteration::Flux_Iterator::
Flux_Iterator(Source_Iteration const &si):
    Vector_Operator(si.phi_size() + si.number_of_augments(),
                    si.phi_size() + si.number_of_augments()),
    si_(si)
{
    
}

void Source_Iteration::Flux_Iterator::
apply(vector<double> &x)
{

    shared_ptr<Ordinate_Sweep_Operator> Linv = dynamic_pointer_cast<Ordinate_Sweep_Operator>(si_.sweeper_);
    shared_ptr<Vector_Operator> D = make_shared<Augmented_Operator>(si_.number_of_augments(), si_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M = make_shared<Augmented_Operator>(si_.number_of_augments(), si_.moment_to_discrete_);
    shared_ptr<Vector_Operator> S = make_shared<Augmented_Operator>(si_.number_of_augments(), si_.scattering_);
    shared_ptr<Vector_Operator> F = make_shared<Augmented_Operator>(si_.number_of_augments(), si_.fission_);

    Linv->include_boundary_source(false);
    
    {
        vector<double> x1(x);
        
        (*S)(x); // moment scattering source
        (*F)(x1); // moment fission source
        
        for (int i = 0; i < si_.phi_size(); ++i)
        {
            x[i] += x1[i];
        }
    }
    
    (*M)(x); // discrete fission+scattering source
    (*Linv)(x); // solution
    (*D)(x); // moment solution
}
