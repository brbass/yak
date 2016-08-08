#ifndef Source_Iteration_hh
#define Source_Iteration_hh

#include <memory>

#include "Solver.hh"
#include "Source_Data.hh"
#include "Vector_Operator.hh"

using std::shared_ptr;

/*
  Uses source iteration to solve the problem
  phi = D Linv M S phi + D Linv q

  D: Discrete to moment
  Linv: Inverse of transport operator
  M: Moment to discrete
  phi: Moment representation of the angular flux
  q: Discrete representation of the internal source
*/
class Source_Iteration : public Solver
{
public:

    // Constructor
    Source_Iteration(int max_iterations,
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
                     shared_ptr<Vector_Operator> fission,
                     shared_ptr<Vector_Operator> preconditioner = shared_ptr<Vector_Operator>());
    
    // Solve fixed source problem
    virtual void solve_steady_state(vector<double> &x) override;

    // Solve k-eigenvalue problem
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x) override;

    // Solve time-dependent problem
    virtual void solve_time_dependent(vector<double> &x) override;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;

    // Size of moment representation of flux
    int phi_size() const
    {
        return moment_to_discrete_->column_size();
    }

    // Number of data values to store for reflection
    int number_of_augments() const
    {
        return source_data_->number_of_augments();
    }

private:

    // Check for convergence based on pointwise error in scalar flux
    bool check_phi_convergence(vector<double> const &x, 
                               vector<double> const &x_old,
                               double &error) const;

    bool check_k_convergence(double k,
                             double k_old,
                             double &error) const;
    
    double update_k(double k_eigenvalue_old,
                    vector<double> const &x,
                    vector<double> const &x_old) const;
    
    bool preconditioned_;
    int max_iterations_;
    int total_iterations_;
    int source_iterations_;
    double tolerance_;

    shared_ptr<Vector_Operator> sweeper_;
    shared_ptr<Vector_Operator> discrete_to_moment_;
    shared_ptr<Vector_Operator> moment_to_discrete_;
    shared_ptr<Vector_Operator> scattering_;
    shared_ptr<Vector_Operator> fission_;
    shared_ptr<Vector_Operator> preconditioner_;
    
    /*
      Computes the first-flight flux, including boundary sources, via source iteration
      b = D Linv q

      D: Discrete to moment
      Linv: Inverse of transport operator
      M: Moment to discrete
      q: Discrete representation of the internal source
      b: Moment representation of the first-flight flux
    */
    class Source_Iterator : public Vector_Operator
    {
    public:
        
        Source_Iterator(Source_Iteration const &source_iteration);
        
    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Source_Iteration const &si_;
    };
    
    /*
      Computes the moments of the flux via source iteration
      phi = D Linv M S phi + b
      
      D: Discrete to moment
      Linv: Inverse of transport operator
      M: Moment to discrete
      phi: Moment representation of the angular flux
      b: Moment representation of the first-flight flux
    */
    class Flux_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Flux_Iterator(Source_Iteration const &source_iteration);
        
    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Source_Iteration const &si_;
    };
};


#endif
