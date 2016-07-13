#ifndef Diffusion_Synthetic_Acceleration_hh
#define Diffusion_Synthetic_Acceleration_hh

#include "Preconditioner.hh"

class AztecOO;
class Epetra_Comm;
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_Operator;
class Epetra_Vector;

class Scattering_Operator;
class Sweep_Operator;

/*
  Performs DSA by either full resolution of the multigroup problem 
  or by assuming that the error in the other groups is zero,
  which allows the use of a removal cross section
*/
class Diffusion_Synthetic_Acceleration : public Preconditioner
{
public:
    
    enum class DSA_Type
    {
        MULTIGROUP,
        REMOVAL
    };

    Diffusion_Synthetic_Acceleration(DSA_Type dsa_type,
                                     shared_ptr<Spatial_Discretization> spatial_discretization,
                                     shared_ptr<Angular_Discretization> angular_discretization,
                                     shared_ptr<Energy_Discretization> energy_discretization,
                                     shared_ptr<Nuclear_Data> nuclear_data,
                                     shared_ptr<Source_Data> source_data,
                                     shared_ptr<Sweep_Operator> sweeper,
                                     shared_ptr<Scattering_Operator> scattering,
                                     shared_ptr<Scattering_Operator> fission);
    
private:
    
    virtual void apply(vector<double> &x) const override;
    virtual void apply_multigroup(vector<double> &x) const;
    virtual void apply_removal(vector<double> &x) const;

    void initialize_trilinos();
    
    int phi_size_;
    int max_iterations_ = 200;
    int kspace_ = 10;
    double tolerance_ = 1e-6;
    DSA_Type dsa_type_;
    shared_ptr<Scattering_Operator> scattering_;
    shared_ptr<Scattering_Operator> fission_;

    shared_ptr<Epetra_Comm> comm_;
    shared_ptr<Epetra_Map> map_;
    shared_ptr<Epetra_Vector> lhs_;
    shared_ptr<Epetra_Vector> rhs_;
    shared_ptr<Epetra_Operator> oper_;
    shared_ptr<Epetra_LinearProblem> problem_;
    shared_ptr<AztecOO> solver_;
};

#endif
