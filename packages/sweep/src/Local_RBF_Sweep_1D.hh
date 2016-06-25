#ifndef Local_RBF_Sweep_1D_hh
#define Local_RBF_Sweep_1D_hh

#include <memory>

#include "Ordinate_Sweep_Operator.hh"

class Amesos_BaseSolver;
class AztecOO;
class Epetra_CrsMatrix;
class Epetra_Comm;
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_Vector;

class Local_RBF_Mesh;

using std::shared_ptr;

/*
  Inverse operator for RBF sweep
  
  Creates and solves the matrix directly
*/
class Local_RBF_Sweep_1D : public Ordinate_Sweep_Operator
{
public:

    // Type of matrix solver
    enum class Solver_Type
    {
        AMESOS,
        AZTECOO
    };

    // Constructor
    Local_RBF_Sweep_1D(shared_ptr<Spatial_Discretization> spatial_discretization,
                       shared_ptr<Angular_Discretization> angular_discretization,
                       shared_ptr<Energy_Discretization> energy_discretization,
                       shared_ptr<Nuclear_Data> nuclear_data,
                       shared_ptr<Source_Data> source_data,
                       Solver_Type solver_type = Solver_Type::AMESOS);
    
protected:

    shared_ptr<Local_RBF_Mesh> rbf_mesh_;

private:
    
    virtual void apply(vector<double> &x) const override;
    void sweep_slab(vector<double> &x) const;
    void initialize_trilinos();
    void add_boundary_point(int b,
                            int i,
                            int o,
                            int g,
                            vector<double> const &x) const;
    void add_internal_point(int i,
                            int o,
                            int g,
                            vector<double> const &x) const;

    int max_iterations_ = 5000;
    double tolerance_ = 1e-6;
    Solver_Type solver_type_;
    
    shared_ptr<Epetra_Comm> comm_;
    shared_ptr<Epetra_Map> map_;
    shared_ptr<Epetra_Vector> lhs_;
    shared_ptr<Epetra_Vector> rhs_;
    shared_ptr<Epetra_CrsMatrix> mat_;
    shared_ptr<Epetra_LinearProblem> problem_;
    shared_ptr<AztecOO> aztec_solver_;
    shared_ptr<Amesos_BaseSolver*> amesos_solver_;
};

#endif
