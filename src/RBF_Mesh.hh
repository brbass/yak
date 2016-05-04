#ifndef RBF_Mesh_hh
#define RBF_Mesh_hh

#include <vector>

#include "Spatial_Discretization.hh"

using std::vector;

class RBF_Mesh : public Spatial_Discretization
{
public:
    
    RBF_Mesh(int dimension,
             int number_of_points,
             Geometry geometry,
             vector<int> const &material,
             vector<double> const &positions);

    virtual int number_of_points()
    {
        return number_of_points_;
    }
    virtual int number_of_cells()
    {
        return number_of_points_;
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
    vector<double> const &point_positions() const
    {
        return point_positions_;
    }

    void check_class_invariants() const;
    virtual void output(pugi::xml_node &output_node);
    
private:
    
    int number_of_points_;

    vector<bool> boundary_nodes_;
    vector<int> boundary_points_;
    vector<int> internal_points_;
    vector<int> material_;
    vector<double> point_positions_;
};

#endif
