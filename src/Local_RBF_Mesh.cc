#include "Local_RBF_Mesh.hh"

Local_RBF_Mesh::
Local_RBF_Mesh(int dimension,
               int number_of_points,
               int number_of_neighbors,
               Geometry geometry,
               Basis_Type basis_type,
               vector<int> const &material,
               vector<double> const &positions,
               vector<double> const &shape_parameter):
    Spatial_Discretization(dimension,
                           geometry),
    RBF_Mesh(dimension,
             number_of_points,
             geometry,
             basis_type,
             material,
             positions,
             shape_parameter)
{
    // Find nearest neighbors
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<int> local_neighbors;
        vector<double> local_distance;
        get_neighbors(i,
                      local_neighbors,
                      local_distance);
        neighbors_.push_back(local_neighbors);
        neighbor_distance_.push_back(local_distance);
    }
    
    // Initialize Trilinos
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        
    }
}

void Local_RBF_Mesh::
convert_to_local(int point,
                 vector<double> &x)
{
    
}

void Local_RBF_Mesh::
get_neighbors(int point,
              vector<int> &local_neighbors,
              vector<int> &local_distance)
{
    local_neighbors.resize(number_of_points_);

    vector<double> distance(number_of_points_);
    shared_ptr<RBF> const local_basis = basis_functions(point);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        shared_ptr<RBF> const other_basis = basis_functions(i);
        vector<double> const other_position = other_basis->position();

        local_neighbors[i] = i;
        distance[i] = local_basis->get_distance_squared(other_position);
    }
    
    partial_sort(local_neighbors.begin(),
                 local_neighbors.begin() + number_of_neighbors_,
                 local_neighbors.end(),
                 [&](int i, int j)
                 {
                     return distance[i] < distance[j];
                 });
    
    local_neighbors.resize(number_of_neigbors_);
    
    local_distance.resize(number_of_neighbors_);

    for (int i = 0; i < number_of_neighbors_; ++i)
    {
        local_distance[i] = distance[local_neighbors[i]];
    }
}
