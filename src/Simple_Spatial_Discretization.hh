#ifndef Simple_Spatial_Discretization_hh
#define Simple_Spatial_Discretization_hh

#include "Spatial_Discretization.hh"

/*
  Simple definition of a spatial discretization for testing
*/
class Simple_Spatial_Discretization : public Spatial_Discretization
{
public:

    Simple_Spatial_Discretization(int dimension,
                                  int number_of_cells,
                                  int number_of_nodes,
                                  int number_of_boundary_cells,
                                  Geometry geometry);

    virtual int number_of_points()
    {
        return number_of_cells_ * number_of_nodes_;
    }

    virtual int number_of_internal_points()
    {
        return (number_of_cells_ - number_of_boundary_cells_) * number_of_nodes_;
    }
    
    virtual int number_of_boundary_points()
    {
        return number_of_boundary_cells_ * number_of_nodes_;
    }

    virtual int number_of_cells()
    {
        return number_of_cells_;
    }

    virtual int number_of_boundary_cells()
    {
        return number_of_boundary_cells_;
    }

    virtual int number_of_nodes()
    {
        return number_of_nodes_;
    }

    virtual vector<bool> const &boundary_nodes() const
    {
        return boundary_nodes_;
    }

    virtual vector<int> const &boundary_cells() const
    {
        return boundary_cells_;
    }

    virtual vector<int> const &internal_cells() const
    {
        return internal_cells_;
    }

    virtual vector<int> const &material() const
    {
        return material_;
    }

    virtual void output(pugi::xml_node &output_node) const
    {
    }

private:

    int number_of_cells_;
    int number_of_nodes_;
    int number_of_boundary_cells_;
    
    vector<bool> boundary_nodes_;
    vector<int> boundary_cells_;
    vector<int> internal_cells_;
    vector<int> material_;
};

#endif
