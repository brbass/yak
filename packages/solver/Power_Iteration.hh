#ifndef Power_Iteration_hh
#define Power_Iteration_hh

#include <memory>

#include "Solver.hh"
#include "Source_Data.hh"
#include "Vector_Operator.hh"

class Epetra_Comm;
class Epetra_Map;

using std::shared_ptr;

/*
  Uses power iteration to solve the problem 
  (I - D Linv M S) phi = D Linv M F phi / k
  
  I: Identity
  D: Discrete to moment
  Linv: Inverse of transport operator
  M: Moment to discrete
  phi: Moment representation of the angular flux
  S: Scattering
  F: Fission
  k: k-eigenvalue
*/
class Power_Iteration : public Solver
{
public:
    
    // Constructor
    Power_Iteration(int max_iterations,
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
    
    // Check for convergence based on pointwise error in scalar flux
    bool check_phi_convergence(vector<double> const &x, 
                               vector<double> const &x_old,
                               double &error);
    
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

    bool check_k_convergence(double k,
                             double k_old,
                             double &error) const;
    
    double update_k(double k_eigenvalue_old,
                    vector<double> const &x,
                    vector<double> const &x_old) const;
    
    bool preconditioned_;
    int max_iterations_;
    int kspace_;
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
      Iteratively computes the moments of the flux using Power methods
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
        Flux_Iterator(Power_Iteration const &power_iteration,
                      bool include_fission = true);
        
    private:
        
        virtual void apply(vector<double> &x) const override;

        bool include_fission_;
        Power_Iteration const &pi_;
    };

    /* 
       Represents the fission operator for an eigenvalue calculation,
       D Linv M F,
    */
    class Fission_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Fission_Iterator(Power_Iteration const &power_iteration);
        
    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Power_Iteration const &pi_;
    };
    
    /* 
       Represents the eigenvalue operator,
       (I - D Linv M S)inv D Linv M F
    */
    class Eigenvalue_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Eigenvalue_Iterator(Power_Iteration const &power_iteration,
                            shared_ptr<Epetra_Comm> comm,
                            shared_ptr<Epetra_Map> map);
        
    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Power_Iteration const &pi_;
        shared_ptr<Epetra_Comm> comm_;
        shared_ptr<Epetra_Map> map_;
    };
};


#endif
