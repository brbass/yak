#include "Simple_Spatial_Discretization.hh"

Simple_Spatial_Discretization::
Simple_Spatial_Discretization(int dimension,
                              int number_of_cells,
                              int number_of_nodes,
                              int number_of_boundary_cells,
                              Geometry geometry):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_cells_(number_of_cells),
    number_of_nodes_(number_of_nodes),
    number_of_materials_(1),
    number_of_boundary_cells_(number_of_boundary_cells)
{
    boundary_nodes_.assign(number_of_boundary_cells * number_of_nodes, true);
    boundary_cells_.resize(number_of_boundary_cells, 0);
    internal_cells_.resize(number_of_cells - number_of_boundary_cells, 0);
    material_.assign(number_of_cells, 0);
    boundary_normal_.assign(number_of_boundary_cells * dimension, 0);
    
    for (int i = 0; i < number_of_boundary_cells; ++i)
    {
        boundary_cells_[i] = i;
        
        boundary_normal_[0 + dimension * i] = 1;
    }
    for (int i = number_of_boundary_cells; i < number_of_cells; ++i)
    {
        internal_cells_[i] = i;
    }
}
