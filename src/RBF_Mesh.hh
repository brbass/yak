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
             vector<double> &positions);

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
    virtual int dimension()
    {
        return dimension_;
    }
    vector<double> &point_positions()
    {
        return point_positions_;
    }

    void check_class_invariants();

private:
    
    int number_of_points_;
    int dimension_;
    
    vector<double> point_positions_;
};

#endif
