#ifndef RBF_Mesh_hh
#define RBF_Mesh_hh

#include <memory>
#include <vector>

#include "RBF.hh"
#include "Spatial_Discretization.hh"

class KD_Tree;
class Solid_Geometry;

using std::shared_ptr;
using std::vector;

/*
  Mesh designed for use with radial basis functions
*/
class RBF_Mesh : public Spatial_Discretization
{
public:

    // Type of basis function
    enum class Basis_Type
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
             int number_of_boundary_points,
             int number_of_internal_points,
             int number_of_neighbors,
             int number_of_materials,
             double shape_multiplier,
             Geometry geometry,
             Basis_Type basis_type,
             vector<int> const &material,
             vector<int> const &boundary_points,
             vector<int> const &internal_points,
             vector<double> const &positions,
             vector<double> const &boundary_normal,
             vector<int> const &surface = vector<int>(),
             vector<int> const &region = vector<int>(),
             shared_ptr<Solid_Geometry> const solid_geometry = shared_ptr<Solid_Geometry>());
    
    RBF_Mesh(int dimension,
             int number_of_points,
             int number_of_boundary_points,
             int number_of_internal_points,
             int number_of_transition_points,
             int number_of_neighbors,
             int number_of_materials,
             double shape_multiplier,
             Geometry geometry,
             Basis_Type basis_type,
             vector<int> const &material,
             vector<int> const &boundary_points,
             vector<int> const &internal_points,
             vector<int> const &transition_points,
             vector<double> const &positions,
             vector<double> const &boundary_normal,
             vector<double> const &transition_normal,
             vector<int> const &surface,
             vector<int> const &region,
             shared_ptr<Solid_Geometry> const solid_geometry);
    
    // Number of points
    virtual int number_of_points() override
    {
        return number_of_points_;
    }

    // Number of points not on the boundary
    virtual int number_of_internal_points() override
    {
        return number_of_internal_points_;
    }
    
    // Number of points on the boundary
    virtual int number_of_boundary_points() override
    {
        return number_of_boundary_points_;
    }

    virtual int number_of_transition_points() override
    {
        return number_of_transition_points_;
    }
    
    // Number of points
    virtual int number_of_cells() override
    {
        return number_of_points_;
    }

    // Number of points on the boundary
    virtual int number_of_boundary_cells() override
    {
        return number_of_boundary_points_;
    }

    // Number of nodes per point: always 1
    virtual int number_of_nodes() override
    {
        return 1;
    }

    // Number of neighbors for each point
    virtual int number_of_neighbors() 
    {
        return number_of_neighbors_;
    }

    virtual int number_of_materials() override
    {
        return number_of_materials_;
    }
    
    virtual Basis_Type basis_type() const
    {
        return basis_type_;
    }
    
    // Which nodes are on the boundary for the boundary cells
    // Should always be true, as number_of_nodes = 1
    virtual vector<bool> const &boundary_nodes() const override
    {
        return boundary_nodes_;
    }

    // Get list of nearest neighbors to a point, sorted by distance
    virtual vector<int> const &neighbors(int i) const
    {
        return neighbors_[i];
    }
    
    // Number of points on the boundary
    virtual vector<int> const &boundary_cells() const override
    {
        return boundary_points_;
    }

    // Number of points not on the boundary
    virtual vector<int> const &internal_cells() const override
    {
        return internal_points_;
    }

    // Number of transition points
    virtual vector<int> const &transition_cells() const override
    {
        return transition_points_;
    }
    
    // Surface the boundary point lies on
    virtual vector<int> const &surface() const
    {
        return surface_;
    }

    // Region the internal point is in
    virtual vector<int> const &region() const
    {
        return region_;
    }

    // Material number for each point
    virtual vector<int> const &material() const override
    {
        return material_;
    }
    
    virtual vector<double> const &boundary_normal() const override
    {
        return boundary_normal_;
    }

    virtual vector<double> const &transition_normal() const override
    {
        return transition_normal_;
    }

    virtual vector<double> const &positions() const
    {
        return point_positions_;
    }
    
    virtual vector<Cell_Type> const &cell_type() const override
    {
        return cell_type_;
    }
    
    // Basis function for a certain point
    virtual shared_ptr<RBF> const &basis_function(int point) const
    {
        return basis_functions_[point];
    }

    virtual shared_ptr<Solid_Geometry> const solid_geometry() const
    {
        return solid_geometry_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const;
    
    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;

    virtual void set_shape_multiplier(double shape_multiplier);
    
protected:

    int number_of_points_;
    int number_of_boundary_points_;
    int number_of_internal_points_;
    int number_of_transition_points_;
    int number_of_neighbors_;
    int number_of_materials_;

    double shape_multiplier_;
    
    Basis_Type basis_type_;
    
    vector<bool> boundary_nodes_;
    vector<int> boundary_points_;
    vector<int> internal_points_;
    vector<int> transition_points_;
    vector<int> material_;
    vector<int> region_;
    vector<int> surface_;
    vector<double> boundary_normal_;
    vector<double> transition_normal_;
    vector<double> point_positions_;
    vector<double> shape_parameter_;
    vector<vector<int> > neighbors_;
    
    vector<Cell_Type> cell_type_;
    
    shared_ptr<KD_Tree> kd_tree_;
    shared_ptr<Solid_Geometry> solid_geometry_;
    
    vector<shared_ptr<RBF> > basis_functions_;
};

#endif
