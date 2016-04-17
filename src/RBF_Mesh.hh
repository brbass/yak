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
    vector<double> const &point_positions() const
    {
        return point_positions_;
    }

    void check_class_invariants() const;
    
private:
    
    int number_of_points_;
    
    vector<double> point_positions_;
};

#endif
