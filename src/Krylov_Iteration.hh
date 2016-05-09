#ifndef Krylov_Iteration_hh
#define Krylov_Iteration_hh

#include <memory>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Solver.hh"
#include "Nuclear_Data.hh"
#include "Source_Data.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

using std::shared_ptr;

class Krylov_Iteration : public Solver
{
public:
    
    Krylov_Iteration(int max_iterations,
                     int kspace,
                     int solver_print,
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
                     shared_ptr<Vector_Operator> fission);
    
    virtual void solve_steady_state(vector<double> &x);
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x);
    virtual void solve_time_dependent(vector<double> &x);
    virtual void output(pugi::xml_node &output_node) const;

    bool check_phi_convergence(vector<double> const &x, 
                               vector<double> const &x_old);
    
    int phi_size() const
    {
        return moment_to_discrete_->column_size();
    }
    int number_of_augments() const
    {
        return source_data_->number_of_augments();
    }

private:
    
    int max_iterations_;
    int kspace_;
    int solver_print_;
    int total_iterations_;
    int source_iterations_;
    double tolerance_;

    shared_ptr<Vector_Operator> sweeper_;
    shared_ptr<Vector_Operator> discrete_to_moment_;
    shared_ptr<Vector_Operator> moment_to_discrete_;
    shared_ptr<Vector_Operator> scattering_;
    shared_ptr<Vector_Operator> fission_;

    class Source_Iterator : public Vector_Operator
    {
    public:
        
        Source_Iterator(Krylov_Iteration const &krylov_iteration);
        
    private:
        
        virtual void apply(vector<double> &x);
        
        Krylov_Iteration const &ki_;
    };
    
    class Flux_Iterator : public Vector_Operator
    {
    public:
        
        Flux_Iterator(Krylov_Iteration const &krylov_iteration);
        
    private:
        
        virtual void apply(vector<double> &x);
        
        Krylov_Iteration const &ki_;
    };
};


#endif
