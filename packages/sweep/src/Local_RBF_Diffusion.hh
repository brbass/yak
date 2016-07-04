#ifndef Local_RBF_Diffusion_hh
#define Local_RBF_Diffusion_hh

#include "Sweep_Operator.hh"

class Local_RBF_Mesh;

class Local_RBF_Diffusion : public Sweep_Operator
{
public:

    // Type of matrix solver
    enum class Solver_Type
    {
        AMESOS,
        AZTECOO
    };

    // Constructor
    Local_RBF_Diffusion(shared_ptr<Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization,
                        shared_ptr<Nuclear_Data> nuclear_data,
                        shared_ptr<Source_Data> source_data,
                        Solver_Type solver_type = Solver_Type::AMESOS);
    
private:

    Solver_Type solver_type_;

    void initialize_trilinos();
    
    virtual void apply(vector<double> &x) const override;

    shared_ptr<Local_RBF_Mesh> rbf_mesh_;
};

#endif
