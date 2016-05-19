#include "Spatial_Discretization_Parser.hh"

using namespace std;

Spatial_Discretization_Parser::
Spatial_Discretization_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node spatial = input_file.child("spatial_discretization");
    
    string spatial_type = child_value<string>(spatial, "type");
    
    if (spatial_type == "finite_element")
    {
        spatial_ = get_fem(spatial);
    }
    else if (spatial_type == "rbf")
    {
        spatial_ = get_rbf(spatial);
    }
    else 
    {
        AssertMsg(false, "spatial discretization type \"" + spatial_type + "\" not found");
    }
}

shared_ptr<Finite_Element_Mesh> Spatial_Discretization_Parser::
get_fem(pugi::xml_node &spatial)
{
    int dimension = child_value<int>(spatial, "dimension");
    int number_of_elements = child_value<int>(spatial, "number_of_elements");
    int number_of_nodes = child_value<int>(spatial, "number_of_nodes");
    string geometry_str = child_value<string>(spatial, "geometry");
    Spatial_Discretization::Geometry geometry;
    
    if (dimension == 1)
    {
        if (geometry_str == "slab")
        {
            geometry = Spatial_Discretization::SLAB;
        }
        else if (geometry_str == "spherical")
        {
            geometry = Spatial_Discretization::SPHERICAL;
        }
        else if (geometry_str == "cylindrical")
        {
            geometry = Spatial_Discretization::CYLINDRICAL;
        }
        else
        {
            AssertMsg(false, "geometry \"" + geometry_str + "\" not found");
        }
    }
    else
    {
        AssertMsg(false, "dimension " + to_string(dimension) + " not available");
    }
    
    pugi::xml_node regions = spatial.child("regions");

    int number_of_regions = child_value<int>(regions, "number_of_regions");
    
    vector<int> material(number_of_elements);
    vector<double> node_positions(number_of_elements * number_of_nodes);

    int cell = 0;
    int x = 0;
    
    for (pugi::xml_node region = regions.child("region"); region; region = region.next_sibling("region"))
    {
        int region_cells = child_value<int>(region, "number_of_elements");
        int region_material = child_value<int>(region, "material");
        double region_length = child_value<double>(region, "length");
        double dx = region_length / region_cells;
        double dn = dx / (number_of_nodes - 1);
        
        for (int i = 0; i < region_cells; ++i)
        {
            material[cell] = region_material;
            
            for (int n = 0; n < number_of_nodes; ++n)
            {
                int k = n + number_of_nodes * cell;

                node_positions[k] = x + dx * i + dn * n;
            }
            
            cell += 1;
        }
        
        x += region_length;
    }
    Assert(cell == number_of_elements);

    string element_type_str = child_value<string>(spatial, "element_type");
    Finite_Element_Mesh::Element_Type element_type;

    if (element_type_str == "cfem")
    {
        element_type = Finite_Element_Mesh::CFEM;
    }
    else
    {
        element_type = Finite_Element_Mesh::DFEM;
    }

    return make_shared<Finite_Element_Mesh>(dimension,
                                            number_of_elements,
                                            number_of_nodes,
                                            geometry,
                                            element_type,
                                            material,
                                            node_positions);
}

shared_ptr<RBF_Mesh> Spatial_Discretization_Parser::
get_rbf(pugi::xml_node &spatial)
{
    int dimension = child_value<int>(spatial, "dimension");
    int number_of_points = child_value<int>(spatial, "number_of_points");
    int local = child_value<int>(spatial, "local");
    double shape_multiplier = child_value<double>(spatial, "shape_multiplier");
    string geometry_str = child_value<string>(spatial, "geometry");
    string basis_str = child_value<string>(spatial, "basis_type");
    Spatial_Discretization::Geometry geometry;
    RBF_Mesh::Basis_Type basis_type;
    
    if (dimension == 1)
    {
        if (geometry_str == "slab")
        {
            geometry = Spatial_Discretization::SLAB;
        }
        else if (geometry_str == "spherical")
        {
            geometry = Spatial_Discretization::SPHERICAL;
        }
        else if (geometry_str == "cylindrical")
        {
            geometry = Spatial_Discretization::CYLINDRICAL;
        }
        else
        {
            AssertMsg(false, "geometry \"" + geometry_str + "\" not found");
        }
    }
    else
    {
        AssertMsg(false, "dimension " + to_string(dimension) + " not available");
    }
    
    if (basis_str == "gaussian")
    {
        basis_type = RBF_Mesh::GAUSSIAN;
    }
    else
    {
        AssertMsg(false, "basis_type \"" + basis_str + "\" not found");
    }

    pugi::xml_node regions = spatial.child("regions");

    int number_of_regions = child_value<int>(regions, "number_of_regions");

    vector<int> material(number_of_points);
    vector<double> positions(number_of_points);
    vector<double> shape_parameter(number_of_points);

    int point = 0;
    int region_x = 0;
    
    for (pugi::xml_node region = regions.child("region"); region; region = region.next_sibling("region"))
    {
        int region_points = child_value<int>(region, "number_of_points");
        int region_material = child_value<int>(region, "material");
        double region_length = child_value<double>(region, "length");
        double dx = region_length / region_points;
        
        if (point == 0)
        {
            material[point] = region_material;
            positions[point] = region_x;
            shape_parameter[point] = shape_multiplier / dx;
            
            point += 1;
        }
        
        for (int i = 0; i < region_points; ++i)
        {
            material[point] = region_material;
            positions[point] = region_x + (0.5 + i) * dx;
            shape_parameter[point] = shape_multiplier / dx;
            
            point += 1;
        }
        
        region_x += region_length;
        
        if (point == number_of_points - 1)
        {
            material[point] = region_material;
            positions[point] = region_x;
            shape_parameter[point] = shape_multiplier / dx;
            
            point += 1;
        }
    }
    Assert(point == number_of_points);

    if (local)
    {
        int number_of_neighbors = child_value<int>(spatial, "number_of_neighbors");
        
        return make_shared<Local_RBF_Mesh>(dimension,
                                           number_of_points,
                                           number_of_neighbors,
                                           geometry,
                                           basis_type,
                                           material,
                                           positions,
                                           shape_parameter);
    }
    else
    {
        return make_shared<RBF_Mesh>(dimension,
                                     number_of_points,
                                     geometry,
                                     basis_type,
                                     material,
                                     positions,
                                     shape_parameter);
    }
}
