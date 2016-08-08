#include "Spatial_Discretization.hh"

Spatial_Discretization::
Spatial_Discretization(int dimension,
                       Geometry geometry):
    dimension_(dimension),
    geometry_(geometry)
{
}

void Spatial_Discretization::
get_cell_type(vector<Cell_Type> &cell_type)
{
    cell_type.resize(number_of_cells());
    
    for (int c = 0; c < number_of_boundary_cells(); ++c)
    {
        int i = boundary_cells()[c];
        
        cell_type[i] = Cell_Type::BOUNDARY;
    }

    int number_of_internal_cells = number_of_cells() - number_of_boundary_cells() - number_of_transition_points();
    
    for (int c = 0; c < number_of_internal_cells; ++c)
    {
        int i = internal_cells()[c];
        
        cell_type[i] = Cell_Type::INTERNAL;
    }
    
    for (int c = 0; c < number_of_transition_points(); ++c)
    {
        int i = transition_cells()[c];

        cell_type[i] = Cell_Type::TRANSITION;
    }
}
