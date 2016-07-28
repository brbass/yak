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

    virtual int number_of_points() override
    {
        return number_of_cells_ * number_of_nodes_;
    }

    virtual int number_of_internal_points() override
    {
        return (number_of_cells_ - number_of_boundary_cells_) * number_of_nodes_;
    }
    
    virtual int number_of_boundary_points() override
    {
        return number_of_boundary_cells_ * number_of_nodes_;
    }

    virtual int number_of_cells() override
    {
        return number_of_cells_;
    }

    virtual int number_of_boundary_cells() override
    {
        return number_of_boundary_cells_;
    }

    virtual int number_of_nodes() override
    {
        return number_of_nodes_;
    }

    virtual int number_of_materials() override
    {
        return number_of_materials_;
    }

    virtual vector<bool> const &boundary_nodes() const override
    {
        return boundary_nodes_;
    }

    virtual vector<int> const &boundary_cells() const override
    {
        return boundary_cells_;
    }

    virtual vector<int> const &internal_cells() const override
    {
        return internal_cells_;
    }

    virtual vector<int> const &material() const override
    {
        return material_;
    }

    virtual vector<double> const &boundary_normal() const override
    {
        return boundary_normal_;
    }
    
    virtual void output(pugi::xml_node &output_node) const override
    {
    }

private:

    int number_of_cells_;
    int number_of_nodes_;
    int number_of_boundary_cells_;
    int number_of_materials_;
    
    vector<bool> boundary_nodes_;
    vector<int> boundary_cells_;
    vector<int> internal_cells_;
    vector<int> material_;
    vector<double> boundary_normal_;
};

#endif
