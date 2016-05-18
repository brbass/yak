#ifndef RBF_Mesh_hh
#define RBF_Mesh_hh

#include <memory>
#include <vector>

#include "RBF.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

class RBF_Mesh : public Spatial_Discretization
{
public:
    
    enum Basis_Type
    {
        GAUSSIAN,
        MULTIQUADRIC,
        INVERSE_MULTIQUADRIC,
        WENDLAND30,
        WENDLAND31,
        WENDLAND32,
        WENDLAND33,
        FALSE
    };

    RBF_Mesh(int dimension,
             int number_of_points,
             Geometry geometry,
             Basis_Type basis_type,
             vector<int> const &material,
             vector<double> const &positions,
             vector<double> const &shape_parameter);

    virtual int number_of_points()
    {
        return number_of_points_;
    }
    virtual int number_of_boundary_points()
    {
        return number_of_boundary_points_;
    }
    virtual int number_of_cells()
    {
        return number_of_points_;
    }
    virtual int number_of_boundary_cells()
    {
        return number_of_boundary_points_;
    }
    virtual int number_of_nodes()
    {
        return 1;
    }
    virtual vector<bool> const &boundary_nodes() const
    {
        return boundary_nodes_;
    }
    virtual vector<int> const &boundary_cells() const
    {
        return boundary_points_;
    }
    virtual vector<int> const &internal_cells() const
    {
        return internal_points_;
    }
    virtual vector<int> const &material() const
    {
        return material_;
    }
    shared_ptr<RBF> const basis_function(int point) const
    {
        return basis_functions_[point];
    }
    
    void check_class_invariants() const;
    virtual void output(pugi::xml_node &output_node) const;
    
protected:
    
    int number_of_points_;
    int number_of_boundary_points_;
    int number_of_internal_points_;

    Basis_Type basis_type_;

    vector<bool> boundary_nodes_;
    vector<int> boundary_points_;
    vector<int> internal_points_;
    vector<int> material_;
    vector<double> point_positions_;
    vector<double> shape_parameter_;
    
    vector<shared_ptr<RBF> > basis_functions_;
};

#endif
