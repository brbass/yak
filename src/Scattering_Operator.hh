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

/*
  Pure virtual class to apply scattering to a moment representation of the flux
*/
class Scattering_Operator : public Vector_Operator
{
public:

    // Types of scattering
    enum Scattering_Type
    {
        COHERENT,
        INCOHERENT,
        FULL
    };

    // Constructor
    Scattering_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization,
                        shared_ptr<Nuclear_Data> nuclear_data,
                        Scattering_Type scattering_type = FULL);
    
protected:

    // Get vector input/output size
    virtual int get_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                         shared_ptr<Angular_Discretization> angular_discretization,
                         shared_ptr<Energy_Discretization> energy_discretization);

    // Type of scattering
    Scattering_Type scattering_type_;

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Nuclear_Data> nuclear_data_;

private: 

    // Apply scattering of chosen type
    virtual void apply(vector<double> &x);

    // Apply within-group and out-of-group scattering
    virtual void apply_full(vector<double> &x) = 0;
    
    // Apply only within-group scattering
    virtual void apply_coherent(vector<double> &x) = 0;
    
    // Apply only out-of-group scattering
    virtual void apply_incoherent(vector<double> &x);
};

#endif
