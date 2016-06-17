#include "Spatial_Discretization_Parser.hh"

using namespace std;

Spatial_Discretization_Parser::
Spatial_Discretization_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node spatial = input_file.child("spatial_discretization");
    
    string spatial_type = XML_Functions::child_value<string>(spatial, "type");
    
    if (spatial_type == "finite_element")
    {
        spatial_ = get_fem(spatial);
    }
    else if (spatial_type == "rbf")
    {
        spatial_ = get_rbf_1d(spatial);
    }
    else if (spatial_type == "rbf_solid")
    {
        spatial_ = get_rbf_solid(spatial);
    }
    else 
    {
        AssertMsg(false, "spatial discretization type \"" + spatial_type + "\" not found");
    }
}

shared_ptr<Finite_Element_Mesh> Spatial_Discretization_Parser::
get_fem(pugi::xml_node &spatial)
{
    int dimension = XML_Functions::child_value<int>(spatial, "dimension");
    int number_of_elements = XML_Functions::child_value<int>(spatial, "number_of_elements");
    int number_of_nodes = XML_Functions::child_value<int>(spatial, "number_of_nodes");
    string geometry_str = XML_Functions::child_value<string>(spatial, "geometry");
    Spatial_Discretization::Geometry geometry;
    
    if (dimension == 1)
    {
        if (geometry_str == "slab")
        {
            geometry = Spatial_Discretization::Geometry::SLAB;
        }
        else if (geometry_str == "spherical")
        {
            geometry = Spatial_Discretization::Geometry::SPHERE;
        }
        else if (geometry_str == "cylindrical")
        {
            geometry = Spatial_Discretization::Geometry::CYLINDER;
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

    int number_of_regions = XML_Functions::child_value<int>(regions, "number_of_regions");
    
    vector<int> material(number_of_elements);
    vector<double> node_positions(number_of_elements * number_of_nodes);

    int cell = 0;
    int x = 0;
    
    for (pugi::xml_node region = regions.child("region"); region; region = region.next_sibling("region"))
    {
        int region_cells = XML_Functions::child_value<int>(region, "number_of_elements");
        int region_material = XML_Functions::child_value<int>(region, "material");
        double region_length = XML_Functions::child_value<double>(region, "length");
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

    string element_type_str = XML_Functions::child_value<string>(spatial, "element_type");
    Finite_Element_Mesh::Element_Type element_type;

    if (element_type_str == "cfem")
    {
        element_type = Finite_Element_Mesh::Element_Type::CFEM;
    }
    else
    {
        element_type = Finite_Element_Mesh::Element_Type::DFEM;
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
get_rbf_1d(pugi::xml_node &spatial)
{
    int dimension = XML_Functions::child_value<int>(spatial, "dimension");
    int number_of_points = XML_Functions::child_value<int>(spatial, "number_of_points");
    int local = XML_Functions::child_value<int>(spatial, "local");
    double shape_multiplier = XML_Functions::child_value<double>(spatial, "shape_multiplier");
    string geometry_str = XML_Functions::child_value<string>(spatial, "geometry");
    string basis_str = XML_Functions::child_value<string>(spatial, "basis_type");
    Spatial_Discretization::Geometry geometry;
    RBF_Mesh::Basis_Type basis_type;
    
    if (dimension == 1)
    {
        if (geometry_str == "slab")
        {
            geometry = Spatial_Discretization::Geometry::SLAB;
        }
        else if (geometry_str == "spherical")
        {
            geometry = Spatial_Discretization::Geometry::SPHERE;
        }
        else if (geometry_str == "cylindrical")
        {
            geometry = Spatial_Discretization::Geometry::CYLINDER;
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
        basis_type = RBF_Mesh::Basis_Type::GAUSSIAN;
    }
    else if (basis_str == "multiquadric")
    {
        basis_type = RBF_Mesh::Basis_Type::MULTIQUADRIC;
    }
    else if (basis_str == "inverse_multiquadric")
    {
        basis_type = RBF_Mesh::Basis_Type::INVERSE_MULTIQUADRIC;
    }
    else if (basis_str == "wendland30")
    {
        basis_type = RBF_Mesh::Basis_Type::WENDLAND30;
    }
    else if (basis_str == "wendland31")
    {
        basis_type = RBF_Mesh::Basis_Type::WENDLAND31;
    }
    else if (basis_str == "wendland32")
    {
        basis_type = RBF_Mesh::Basis_Type::WENDLAND32;
    }
    else if (basis_str == "wendland33")
    {
        basis_type = RBF_Mesh::Basis_Type::WENDLAND33;
    }
    else
    {
        AssertMsg(false, "basis_type \"" + basis_str + "\" not found");
    }

    pugi::xml_node regions = spatial.child("regions");

    int number_of_regions = XML_Functions::child_value<int>(regions, "number_of_regions");

    vector<int> material(number_of_points);
    vector<double> positions(number_of_points);
    vector<double> shape_parameter(number_of_points);

    int point = 0;
    int region_x = 0;
    
    for (pugi::xml_node region = regions.child("region"); region; region = region.next_sibling("region"))
    {
        int region_points = XML_Functions::child_value<int>(region, "number_of_points");
        int region_material = XML_Functions::child_value<int>(region, "material");
        double region_length = XML_Functions::child_value<double>(region, "length");
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

    int number_of_boundary_points = 2;
    int number_of_internal_points = number_of_points - 2;
    vector<int> boundary_points;
    boundary_points.push_back(0);
    boundary_points.push_back(number_of_points - 1);
    
    vector<int> internal_points(number_of_internal_points);
    for (int i = 1; i < number_of_points - 1; ++i)
    {
        internal_points[i - 1].push_back(i);
    }
    vector<double> boundary_normal = {-1, 1};
    
    if (local)
    {
        int number_of_neighbors = XML_Functions::child_value<int>(spatial, "number_of_neighbors");
        
        return make_shared<Local_RBF_Mesh>(dimension,
                                           number_of_points,
                                           number_of_boundary_points,
                                           number_of_internal_points,
                                           number_of_neighbors,
                                           geometry,
                                           basis_type,
                                           material,
                                           boundary_points,
                                           internal_points,
                                           positions,
                                           shape_parameter,
                                           boundary_normal);
    }
    else
    {
        return make_shared<RBF_Mesh>(dimension,
                                     number_of_points,
                                     number_of_boundary_points,
                                     number_of_internal_points,
                                     geometry,
                                     basis_type,
                                     material,
                                     boundary_points,
                                     internal_points,
                                     positions,
                                     shape_parameter,
                                     boundary_normal);
    }
}

shared_ptr<RBF_Mesh> Spatial_Discretization_Parser::
get_rbf_solid(pugi::xml_node &spatial)
{
    Solid_Geometry_Parser solid_geometry_parser(spatial);

    shared_ptr<Solid_Geometry> solid_geometry = solid_geometry_parser.get_ptr();
    
    
}

void Spatial_Discretization_Parser::
get_solid_points(shared_ptr<Solid_Geometry> solid_geometry,
                 pugi::xml_node &spatial,
                 int &number_of_points,
                 int &number_of_boundary_points,
                 int &number_of_internal_points,
                 vector<int> &material,
                 vector<int> &boundary_points,
                 vector<int> &internal_points,
                 vector<double> &positions,
                 vector<double> &boundary_normal)
{
    int dimension = solid_geometry->dimension();
    int number_of_points = ;
    int number_of_boundary_points = ;
    int number_of_internal_points = ;
    
    material.resize(number_of_points);
    boundary_points.resize(number_of_boundary_points);
    internal_points.resize(number_of_internal_points);
    positions.resize(number_of_points * dimension);
    boundary_normal.resize(number_of_boundary_points * dimension);
    
    // Get parameters for bounding sphere
    
    double bounding_radius = child_value<double>(spatial, "bounding_radius");
    vector<double> bounding_origin = child_value<double>(spatial, "bounding_origin");

    int current_boundary_point = 0;
    while (current_boundary_point < number_of_boundary_points)
    {
        
    }
    switch(dimension)
    {
    case 2:
    {
        
    }
    case 3:
    {
        
    }
    default:
        AssertMsg(false, "dimension not found");
    }
}
