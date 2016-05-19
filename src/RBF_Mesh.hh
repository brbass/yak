#ifndef RBF_Mesh_hh
#define RBF_Mesh_hh

#include <memory>
#include <vector>

#include "RBF.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

/*
  Mesh designed for use with radial basis functions
*/
class RBF_Mesh : public Spatial_Discretization
{
public:

    // Type of basis function
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

    // Constructor
    RBF_Mesh(int dimension,
             int number_of_points,
             Geometry geometry,
             Basis_Type basis_type,
             vector<int> const &material,
             vector<double> const &positions,
             vector<double> const &shape_parameter);

    // Number of points
    virtual int number_of_points()
    {
        return number_of_points_;
    }

    // Number of points on the boundary
    virtual int number_of_boundary_points()
    {
        return number_of_boundary_points_;
    }

    // Number of points
    virtual int number_of_cells()
    {
        return number_of_points_;
    }

    // Number of points on the boundary
    virtual int number_of_boundary_cells()
    {
        return number_of_boundary_points_;
    }

    // Number of nodes per point: always 1
    virtual int number_of_nodes()
    {
        return 1;
    }

    // Which nodes are on the boundary for the boundary cells
    // Should always be true, as number_of_nodes = 1
    virtual vector<bool> const &boundary_nodes() const
    {
        return boundary_nodes_;
    }

    // Number of points on the boundary
    virtual vector<int> const &boundary_cells() const
    {
        return boundary_points_;
    }

    // Number of points not on the boundary
    virtual vector<int> const &internal_cells() const
    {
        return internal_points_;
    }

    // Material number for each point
    virtual vector<int> const &material() const
    {
        return material_;
    }

    virtual vector<double> const &boundary_normal() const
    {
        return boundary_normal_;
    }
    
    // Basis function for a certain point
    shared_ptr<RBF> const basis_function(int point) const
    {
        return basis_functions_[point];
    }

    // Check class invariants
    void check_class_invariants() const;

    // Output data to XML file
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
    vector<double> boundary_normal_;
    vector<double> point_positions_;
    vector<double> shape_parameter_;
    
    vector<shared_ptr<RBF> > basis_functions_;
};

#endif
