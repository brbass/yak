#ifndef Sweep_Operator_hh
#define Sweep_Operator_hh

#include <memory>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Nuclear_Data;
class Source_Data;
class Spatial_Discretization;

using std::shared_ptr;

/*
  Transport inverse operator
  
  \Omega \cdot \nabla \psi + \Sigma_t \psi = q
*/
class Sweep_Operator : public Vector_Operator
{
public:

    enum class Sweep_Type
    {
        MOMENT,
        ORDINATE
    };
    
    // Constructor
    Sweep_Operator(Sweep_Type sweep_type,
                   shared_ptr<Spatial_Discretization> spatial_discretization,
                   shared_ptr<Angular_Discretization> angular_discretization,
                   shared_ptr<Energy_Discretization> energy_discretization,
                   shared_ptr<Nuclear_Data> nuclear_data,
                   shared_ptr<Source_Data> source_data);
    
    // Include boundary source in sweep
    virtual void include_boundary_source(bool include_source)
    {
        include_boundary_source_ = include_source;
    }

    // Sweep type
    virtual Sweep_Type sweep_type() const
    {
        return sweep_type_;
    }

    virtual void use_removal(bool use)
    {
        use_removal_ = use;
    }
    
protected:

    bool include_boundary_source_;
    bool use_removal_;
    Sweep_Type sweep_type_;
    
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
    shared_ptr<Source_Data> source_data_;

private:
    
    virtual void apply(vector<double> &x) const override = 0;
};

#endif
