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

/*
  Uses Krylov methods to solve the problem 
  (I - D Linv M S) phi = D Linv q
  
  I: Identity
  D: Discrete to moment
  Linv: Inverse of transport operator
  M: Moment to discrete
  phi: Moment representation of the angular flux
  q: Discrete representation of the internal source
*/
class Krylov_Iteration : public Solver
{
public:

    // Constructor
    Krylov_Iteration(int max_iterations,
                     int kspace, // Number of past guesses to store
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

    // Solve fixed source problem
    virtual void solve_steady_state(vector<double> &x);

    // Solve k-eigenvalue problem
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x);

    // Solve time-dependent problem
    virtual void solve_time_dependent(vector<double> &x);

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const;

    // Check for convergence based on pointwise error in scalar flux
    bool check_phi_convergence(vector<double> const &x, 
                               vector<double> const &x_old);

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

        // Constructor
        Source_Iterator(Krylov_Iteration const &krylov_iteration);
        
    private:
        
        virtual void apply(vector<double> &x);
        
        Krylov_Iteration const &ki_;
    };

    /*
      Iteratively computes the moments of the flux using Krylov methods
      (I - D Linv M S) phi = b

      I: Identity
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
        Flux_Iterator(Krylov_Iteration const &krylov_iteration);
        
    private:
        
        virtual void apply(vector<double> &x);
        
        Krylov_Iteration const &ki_;
    };
};


#endif
