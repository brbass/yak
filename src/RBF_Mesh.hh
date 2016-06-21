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
             Geometry geometry,
             Basis_Type basis_type,
             vector<int> const &material,
             vector<int> const &boundary_points,
             vector<int> const &internal_points,
             vector<double> const &positions,
             vector<double> const &shape_parameter,
             vector<double> const &boundary_normal);

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

    // Which nodes are on the boundary for the boundary cells
    // Should always be true, as number_of_nodes = 1
    virtual vector<bool> const &boundary_nodes() const override
    {
        return boundary_nodes_;
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

    // Material number for each point
    virtual vector<int> const &material() const override
    {
        return material_;
    }

    virtual vector<double> const &boundary_normal() const override
    {
        return boundary_normal_;
    }

    // Basis function for a certain point
    virtual shared_ptr<RBF> const &basis_function(int point) const
    {
        return basis_functions_[point];
    }

    // Check class invariants
    virtual void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;

    // virtual void get_neighbors(int point,
    //                            vector<int> &local_neighbors);
    
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

    // shared_ptr<KD_Adaptor> kd_adaptor_;
    // shared_ptr<KDTreeSingleIndexAdaptor<L2_Adaptor<double, KD_Adaptor> > > 
    
    // class KD_Adaptor
    // {
    // public:
        
    //     KD_Adaptor(RBF_Mesh const &rbf_mesh);

    //     inline int kdtree_get_point_count() const
    //     {
    //         return rbf_mesh_.number_of_points();
    //     }

    //     inline double kdtree_get_pt(const int idx, int dim) const
    //     {
    //         return rbf_mesh_.basis_functions_[idx]->position()[dim];
    //     }

    //     inline double kdtree_distance(const double *p1, const int idx_p2, int size)
    //     {
    //         vector<double> const p1_vec(p1, p1 + size);
        
    //         return rbf_mesh_.basis_functions_[idx]->get_distance_squared(p1_vec);
    //     }
        
    //     template <class BBOX>
    //     bool kdtree_get_bbox(BBOX &/*bb*/) const
    //     {
    //         return false;
    //     }
        
    // private:
        
    //     RBF_Mesh const &rbf_mesh_;
    // }
};

#endif
