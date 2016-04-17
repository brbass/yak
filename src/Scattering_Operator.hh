#ifndef Scattering_Operator_hh
#define Scattering_Operator_hh

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Nuclear_Data.hh"
#include "Spatial_Discretization.hh"
#include "Vector_Operator.hh"

using std::shared_ptr;
using std::vector;

class Scattering_Operator : public Vector_Operator
{
public:
    
    enum Scattering_Type
    {
        COHERENT,
        INCOHERENT,
        FULL
    }

    Scattering_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization,
                        shared_ptr<Nuclear_Data> nuclear_data,
                        scattering_type = FULL);
    
private: 

    virtual void apply(vector<double> &x);
    
    virtual void apply_full(vector<double> &x) = 0;
    virtual void apply_coherent(vector<double> &x) = 0;
    virtual void apply_incoherent(vector<double> &x);
    
    Scattering_Type scattering_type_;
    
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;
}

#endif
