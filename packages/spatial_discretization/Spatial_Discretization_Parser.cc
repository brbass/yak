#include "Spatial_Discretization_Parser.hh"

#include "Finite_Element_Mesh.hh"
#include "Local_RBF_Mesh.hh"
#include "Random_Number_Generator.hh"
#include "RBF_Mesh.hh"
#include "Solid_Geometry.hh"
#include "Solid_Geometry_Parser.hh"
#include "Vector_Functions_2D.hh"
#include "Vector_Functions_3D.hh"

using namespace std;
namespace vf2 = Vector_Functions_2D;
namespace vf3 = Vector_Functions_3D;

namespace // anonymous
{
    Random_Number_Generator<double> mu_rng(-1, // lower
                                           1, // upper
                                           0); // seed
    Random_Number_Generator<double> theta_rng(0, // lower
                                              2 * M_PI, // upper
                                              1); // seed
    
    // Get a random point inside a cube
    
    void get_point(int dimension,
                   double bounding_radius,
                   vector<double> const &bounding_origin,
                   vector<double> &point)
    {
        point.resize(dimension);
        
        for (int d = 0; d < dimension; ++d)
        {
            point[d] = bounding_radius * mu_rng.scalar() + bounding_origin[d];
        }
        
        // switch(dimension)
        // {
        // case 2:
        // {
        //     double theta = 2 * M_PI * rng.random_double();
        //     double r = bounding_radius * rng.random_double();

        //     point[0] = bounding_origin[0] + r * cos(theta);
        //     point[1] = bounding_origin[1] + r * sin(theta);
            
        //     break;
        // }
        // case 3:
        // {
        //     double theta = 2 * M_PI * rng.random_double();
        //     double r = bounding_radius * rng.random_double();
        //     double mu = 2 * rng.random_double() - 1
        //     double sqrt_mu = sqrt(1 - mu * mu);
            
        //     point[0] = bounding_origin[0] + r * sqrt_mu * cos(theta);
        //     point[1] = bounding_origin[1] + r * sqrt_mu * sin(theta);
        //     point[2] = bounding_origin[2] + r * mu;

        //     break;
        // }
        // }
    }

    // Get a random, inward-facing ray from the edge of a sphere
    void get_ray(int dimension,
                 double bounding_radius,
                 vector<double> const &bounding_origin,
                 vector<double> &origin,
                 vector<double> &direction)
    {
        double mu = mu_rng.scalar();
        double theta = theta_rng.scalar();

        vector<double> normal(dimension, 0);
        
        origin.assign(dimension, 0);
        direction.assign(dimension, 0);
        
        switch(dimension)
        {
        case 2:
        {
            origin[0] = bounding_origin[0] + bounding_radius * cos(theta);
            origin[1] = bounding_origin[1] + bounding_radius * sin(theta);
            
            normal[0] = cos(theta);
            normal[1] = sin(theta);

            direction[0] = normal[0];
            direction[1] = normal[1];

            while (vf2::dot(direction, normal) > 0)
            {
                double phi = theta_rng.scalar();
                
                direction[0] = cos(phi);
                direction[1] = sin(phi);
            }
            
            break;
        }
        case 3:
        {
            double sqrt_mu = sqrt(1 - mu * mu);
            
            origin[0] = bounding_origin[0] + bounding_radius * sqrt_mu * cos(theta);
            origin[1] = bounding_origin[1] + bounding_radius * sqrt_mu * sin(theta);
            origin[2] = bounding_origin[2] + bounding_radius * mu;

            normal[0] = sqrt_mu * cos(theta);
            normal[1] = sqrt_mu * sin(theta);
            normal[2] = mu;

            direction[0] = normal[0];
            direction[1] = normal[1];
            direction[2] = normal[2];

            while (vf3::dot(direction, normal) > 0)
            {
                double xi = mu_rng.scalar();
                double phi = theta_rng.scalar();
                double sqrt_xi = sqrt(1 - xi * xi);
                
                direction[0] = sqrt_xi * cos(phi);
                direction[1] = sqrt_xi * sin(phi);
                direction[2] = xi;
            }
            
            break;
        }
        default:
            AssertMsg(false, "dimension not found");
        }
        
    }

    // Check whether point is too close to others
    
    bool point_works(int num_points,
                     int dimension,
                     double min_distance,
                     vector<double> const &point,
                     vector<double> const &positions)
    {
        for (int i = 0; i < num_points; ++i)
        {
            double distance = 0;

            for (int d = 0; d < dimension; ++d)
            {
                double k = point[d] - positions[d + dimension * i];
                
                distance += k * k;
            }

            distance = sqrt(distance);

            if (distance < min_distance)
            {
                return false;
            }
        }

        return true;
    }
    
    void get_solid_points(shared_ptr<Solid_Geometry> solid_geometry,
                          pugi::xml_node &spatial,
                          int &dimension,
                          int &number_of_points,
                          int &number_of_boundary_points,
                          int &number_of_internal_points,
                          vector<int> &surfaces,
                          vector<int> &regions,
                          vector<int> &material,
                          vector<int> &boundary_points,
                          vector<int> &internal_points,
                          vector<double> &positions,
                          vector<double> &boundary_normal)
    {
        dimension = solid_geometry->dimension();

        surfaces.resize(0);
        regions.resize(0);
        material.resize(0);
        boundary_points.resize(0);
        internal_points.resize(0);
        positions.resize(0);
        boundary_normal.resize(0);
        
        // Get parameters for bounding sphere
        
        int max_attempts = XML_Functions::child_value<int>(spatial, "max_attempts");
        double min_distance_boundary = XML_Functions::child_value<double>(spatial, "min_distance_boundary");
        double min_distance_internal = XML_Functions::child_value<double>(spatial, "min_distance_internal");
        double bounding_radius = XML_Functions::child_value<double>(spatial, "bounding_radius");
        vector<double> bounding_origin = XML_Functions::child_vector<double>(spatial, "bounding_origin", dimension);
        
        // Find boundary points

        int current_boundary_point = 0;
        int current_point = 0;
        int current_internal_point = 0;
        int num_attempts = 0;
        
        while (num_attempts < max_attempts)
        {
            vector<double> position;
            vector<double> direction;
            vector<double> normal;
            
            int surface = Solid_Geometry::NO_SURFACE;
            int boundary_region = Solid_Geometry::NO_REGION;
            double distance;
            
            // Find ray that intersects with problem
            
            while (surface == Solid_Geometry::NO_SURFACE)
            {
                // Get random ray in random inward-facing direction
                
                get_ray(dimension,
                        bounding_radius,
                        bounding_origin,
                        position,
                        direction);
                
                // See if ray intersects with a problem boundary
                // Move particle to the problem boundary
                
                surface = solid_geometry->next_boundary(position,
                                                        direction,
                                                        boundary_region,
                                                        distance,
                                                        position);
            }
            
            // Find all intersections of the ray
            
            while(surface != Solid_Geometry::NO_SURFACE)
            {
                // Check whether point is too close to others
                
                if (point_works(current_point,
                                dimension,
                                min_distance_boundary,
                                position,
                                positions))
                {
                    surfaces.push_back(surface);
                    Check(solid_geometry->surface(surface)->normal_direction(position,
                                                                             normal));
                    
                    material.push_back(solid_geometry->region(boundary_region)->material());
                    boundary_points.push_back(current_point);
                    
                    for (int d = 0; d < dimension; ++d)
                    {
                        positions.push_back(position[d]);
                        boundary_normal.push_back(normal[d]);
                    }
                    
                    current_boundary_point += 1;
                    current_point += 1;
                }
                else
                {
                    num_attempts += 1;
                }
                
                // Move particle ahead slightly to avoid finding the same boundary again
                
                solid_geometry->new_position(solid_geometry->delta_distance(),
                                             position,
                                             direction,
                                             position);
                
                // Find the next boundary
                
                surface = solid_geometry->next_boundary(position,
                                                        direction,
                                                        boundary_region,
                                                        distance,
                                                        position);
            }
        }
        
        // Get internal points
        
        num_attempts = 0;
        while (num_attempts < max_attempts)
        {
            int region = Solid_Geometry::NO_REGION;
            vector<double> point;
            
            while(region == Solid_Geometry::NO_REGION)
            {
                // Find random point

                get_point(dimension,
                          bounding_radius,
                          bounding_origin,
                          point);
                
                region = solid_geometry->find_region(point);
            }
            
            if (point_works(current_point,
                            dimension,
                            min_distance_internal,
                            point,
                            positions))
            {
                regions.push_back(region);
                material.push_back(solid_geometry->region(region)->material());
                internal_points.push_back(current_point);

                for (int d = 0; d < dimension; ++d)
                {
                    positions.push_back(point[d]);
                }
                
                current_point += 1;
                current_internal_point += 1;
            }
            else
            {
                num_attempts += 1;
            }
        }
        
        number_of_points = current_point;
        number_of_boundary_points = current_boundary_point;
        number_of_internal_points = current_internal_point;
    }

    void get_solid_points(shared_ptr<Solid_Geometry> solid_geometry,
                          pugi::xml_node &spatial,
                          int &dimension,
                          int &number_of_materials,
                          int &number_of_points,
                          int &number_of_boundary_points,
                          int &number_of_internal_points,
                          int &number_of_transition_points,
                          vector<int> &surfaces,
                          vector<int> &regions,
                          vector<int> &material,
                          vector<int> &boundary_points,
                          vector<int> &internal_points,
                          vector<int> &transition_points,
                          vector<double> &positions,
                          vector<double> &boundary_normal,
                          vector<double> &transition_normal)
    {
        dimension = solid_geometry->dimension();
        number_of_materials = solid_geometry->number_of_materials();
        
        surfaces.resize(0);
        regions.resize(0);
        material.resize(0);
        boundary_points.resize(0);
        internal_points.resize(0);
        transition_points.resize(0);
        positions.resize(0);
        boundary_normal.resize(0);
        transition_normal.resize(0);
        
        // Get parameters for bounding sphere
        
        int max_attempts = XML_Functions::child_value<int>(spatial, "max_attempts");
        double min_distance_boundary = XML_Functions::child_value<double>(spatial, "min_distance_boundary");
        double min_distance_internal = XML_Functions::child_value<double>(spatial, "min_distance_internal");
        double bounding_radius = XML_Functions::child_value<double>(spatial, "bounding_radius");
        vector<double> bounding_origin = XML_Functions::child_vector<double>(spatial, "bounding_origin", dimension);
        
        // Find boundary points

        int current_boundary_point = 0;
        int current_point = 0;
        int current_internal_point = 0;
        int current_transition_point = 0;
        int num_attempts = 0;
        
        while (num_attempts < max_attempts)
        {
            vector<double> position;
            vector<double> direction;
            vector<double> normal;
            
            int surface = Solid_Geometry::NO_SURFACE;
            int boundary_region = Solid_Geometry::NO_REGION;
            double distance;
            
            // Find ray that intersects with problem
            
            while (surface == Solid_Geometry::NO_SURFACE)
            {
                // Get random ray in random inward-facing direction
                
                get_ray(dimension,
                        bounding_radius,
                        bounding_origin,
                        position,
                        direction);
                
                // See if ray intersects with a problem boundary
                // Move particle to the problem boundary
                
                surface = solid_geometry->next_intersection(position,
                                                            direction,
                                                            boundary_region,
                                                            distance,
                                                            position);
            }
            
            // Find all intersections of the ray
            
            while(surface != Solid_Geometry::NO_SURFACE)
            {
                // Check whether point is too close to others
                
                if (point_works(current_point,
                                dimension,
                                min_distance_boundary,
                                position,
                                positions))
                {
                    surfaces.push_back(surface);
                    Check(solid_geometry->surface(surface)->normal_direction(position,
                                                                             normal));

                    vector<double> position_positive;
                    vector<double> position_negative;
                    double delta_distance = solid_geometry->delta_distance();
                        
                    solid_geometry->new_position(delta_distance,
                                                 position,
                                                 normal,
                                                 position_positive);
                    solid_geometry->new_position(-delta_distance,
                                                 position,
                                                 normal,
                                                 position_negative);
                    
                    int region_positive = solid_geometry->find_region(position_positive);
                    int region_negative = solid_geometry->find_region(position_negative);

                    if (region_positive == Solid_Geometry::Geometry_Errors::NO_REGION)
                    {
                        material.push_back(solid_geometry->region(region_negative)->material());
                        boundary_points.push_back(current_point);
                            
                        for (int d = 0; d < dimension; ++d)
                        {
                            positions.push_back(position[d]);
                            boundary_normal.push_back(normal[d]);
                        }
                        
                        current_boundary_point += 1;
                    }
                    else if (region_negative == Solid_Geometry::Geometry_Errors::NO_REGION)
                    {
                        material.push_back(solid_geometry->region(region_positive)->material());
                        boundary_points.push_back(current_point);
                        
                        for (int d = 0; d < dimension; ++d)
                        {
                            positions.push_back(position[d]);
                            boundary_normal.push_back(normal[d]);
                        }
                        
                        current_boundary_point += 1;
                    }
                    else
                    {
                        int material_positive = solid_geometry->region(region_positive)->material();
                        int material_negative = solid_geometry->region(region_negative)->material();
                        material.push_back(material_positive + number_of_materials * material_negative);
                        transition_points.push_back(current_point);
                        
                        for (int d = 0; d < dimension; ++d)
                        {
                            positions.push_back(position[d]);
                            transition_normal.push_back(normal[d]);
                        }
                        
                        current_transition_point += 1;
                    }
                    
                    current_point += 1;
                }
                else
                {
                    num_attempts += 1;
                }
                
                // Move particle ahead slightly to avoid finding the same boundary again
                
                solid_geometry->new_position(solid_geometry->delta_distance(),
                                             position,
                                             direction,
                                             position);
                
                // Find the next boundary
                
                surface = solid_geometry->next_intersection(position,
                                                            direction,
                                                            boundary_region,
                                                            distance,
                                                            position);
            }
        }
        
        // Get internal points
        
        num_attempts = 0;
        while (num_attempts < max_attempts)
        {
            int region = Solid_Geometry::NO_REGION;
            vector<double> point;
            
            while(region == Solid_Geometry::NO_REGION)
            {
                // Find random point

                get_point(dimension,
                          bounding_radius,
                          bounding_origin,
                          point);
                
                region = solid_geometry->find_region(point);
            }
            
            if (point_works(current_point,
                            dimension,
                            min_distance_internal,
                            point,
                            positions))
            {
                regions.push_back(region);
                material.push_back(solid_geometry->region(region)->material());
                internal_points.push_back(current_point);

                for (int d = 0; d < dimension; ++d)
                {
                    positions.push_back(point[d]);
                }
                
                current_point += 1;
                current_internal_point += 1;
            }
            else
            {
                num_attempts += 1;
            }
        }
        
        number_of_points = current_point;
        number_of_transition_points = current_transition_point;
        number_of_boundary_points = current_boundary_point;
        number_of_internal_points = current_internal_point;
    }
} // namespace

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
            
            point += 1;
        }
        
        for (int i = 0; i < region_points; ++i)
        {
            material[point] = region_material;
            positions[point] = region_x + (0.5 + i) * dx;
            
            point += 1;
        }
        
        region_x += region_length;
        
        if (point == number_of_points - 1)
        {
            material[point] = region_material;
            positions[point] = region_x;
            
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
        internal_points[i - 1] = i;
    }
    vector<double> boundary_normal = {-1, 1};
    
    int number_of_neighbors = XML_Functions::child_value<int>(spatial, "number_of_neighbors");
        
    if (local)
    { 
        string coefficient_str = XML_Functions::child_value<string>(spatial, "coefficient_type");
        Local_RBF_Mesh::Coefficient_Type coefficient_type;
        
        if (coefficient_str == "alpha")
        {
            coefficient_type = Local_RBF_Mesh::Coefficient_Type::ALPHA;
        }
        else if (coefficient_str == "phi")
        {
            coefficient_type = Local_RBF_Mesh::Coefficient_Type::PHI;
        }
        else
        {
            AssertMsg(false, "coefficient type \"" + coefficient_str + "\" not found");
        }
        
        return make_shared<Local_RBF_Mesh>(dimension,
                                           number_of_points,
                                           number_of_boundary_points,
                                           number_of_internal_points,
                                           number_of_neighbors,
                                           shape_multiplier,
                                           geometry,
                                           basis_type,
                                           coefficient_type,
                                           material,
                                           boundary_points,
                                           internal_points,
                                           positions,
                                           boundary_normal);
    }
    else
    {
        return make_shared<RBF_Mesh>(dimension,
                                     number_of_points,
                                     number_of_boundary_points,
                                     number_of_internal_points,
                                     number_of_neighbors,
                                     shape_multiplier,
                                     geometry,
                                     basis_type,
                                     material,
                                     boundary_points,
                                     internal_points,
                                     positions,
                                     boundary_normal);
    }
}

shared_ptr<RBF_Mesh> Spatial_Discretization_Parser::
get_rbf_solid(pugi::xml_node &spatial)
{
    Solid_Geometry_Parser solid_geometry_parser(spatial);

    shared_ptr<Solid_Geometry> solid_geometry = solid_geometry_parser.get_ptr();
    
    int dimension;
    int number_of_points;
    int number_of_materials = solid_geometry->number_of_materials();
    int number_of_boundary_points;
    int number_of_internal_points;
    int number_of_transition_points = 0;
    vector<int> material;
    vector<int> boundary_points;
    vector<int> internal_points;
    vector<int> transition_points;
    vector<int> surface;
    vector<int> region;
    vector<double> positions;
    vector<double> boundary_normal;
    vector<double> transition_normal;

    int include_transition_points = XML_Functions::child_value<int>(spatial, "include_transition_points");
    
    if (include_transition_points)
    {
        get_solid_points(solid_geometry,
                         spatial,
                         dimension,
                         number_of_materials,
                         number_of_points,
                         number_of_boundary_points,
                         number_of_internal_points,
                         number_of_transition_points,
                         surface,
                         region,
                         material,
                         boundary_points,
                         internal_points,
                         transition_points,
                         positions,
                         boundary_normal,
                         transition_normal);
    }
    else
    {
        get_solid_points(solid_geometry,
                         spatial,
                         dimension,
                         number_of_points,
                         number_of_boundary_points,
                         number_of_internal_points,
                         surface,
                         region,
                         material,
                         boundary_points,
                         internal_points,
                         positions,
                         boundary_normal);
    }
    
    int local = XML_Functions::child_value<int>(spatial, "local");
    int number_of_neighbors = XML_Functions::child_value<int>(spatial, "number_of_neighbors");
    double shape_multiplier = XML_Functions::child_value<double>(spatial, "shape_multiplier");
    string basis_str = XML_Functions::child_value<string>(spatial, "basis_type");
    
    Spatial_Discretization::Geometry geometry = Spatial_Discretization::Geometry::CARTESIAN;
    
    RBF_Mesh::Basis_Type basis_type;
    
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

    if (local)
    {
        string coefficient_str = XML_Functions::child_value<string>(spatial, "coefficient_type");
        Local_RBF_Mesh::Coefficient_Type coefficient_type;
        
        if (coefficient_str == "alpha")
        {
            coefficient_type = Local_RBF_Mesh::Coefficient_Type::ALPHA;
        }
        else if (coefficient_str == "phi")
        {
            coefficient_type = Local_RBF_Mesh::Coefficient_Type::PHI;
        }
        else
        {
            AssertMsg(false, "coefficient type \"" + coefficient_str + "\" not found");
        }

        if (include_transition_points)
        {
            return make_shared<Local_RBF_Mesh>(dimension,
                                               number_of_points,
                                               number_of_boundary_points,
                                               number_of_internal_points,
                                               number_of_transition_points,
                                               number_of_neighbors,
                                               number_of_materials,
                                               shape_multiplier,
                                               geometry,
                                               basis_type,
                                               coefficient_type,
                                               material,
                                               boundary_points,
                                               internal_points,
                                               transition_points,
                                               positions,
                                               boundary_normal,
                                               transition_normal,
                                               surface,
                                               region,
                                               solid_geometry);
        }
        else
        {
            return make_shared<Local_RBF_Mesh>(dimension,
                                               number_of_points,
                                               number_of_boundary_points,
                                               number_of_internal_points,
                                               number_of_neighbors,
                                               shape_multiplier,
                                               geometry,
                                               basis_type,
                                               coefficient_type,
                                               material,
                                               boundary_points,
                                               internal_points,
                                               positions,
                                               boundary_normal,
                                               surface,
                                               region,
                                               solid_geometry);
        }
    }
    else
    {
        if (include_transition_points)
        {
            return make_shared<RBF_Mesh>(dimension,
                                         number_of_points,
                                         number_of_boundary_points,
                                         number_of_internal_points,
                                         number_of_transition_points,
                                         number_of_neighbors,
                                         number_of_materials,
                                         shape_multiplier,
                                         geometry,
                                         basis_type,
                                         material,
                                         boundary_points,
                                         internal_points,
                                         transition_points,
                                         positions,
                                         boundary_normal,
                                         transition_normal,
                                         surface,
                                         region,
                                         solid_geometry);
        }
        else
        {
            return make_shared<RBF_Mesh>(dimension,
                                         number_of_points,
                                         number_of_boundary_points,
                                         number_of_internal_points,
                                         number_of_neighbors,
                                         shape_multiplier,
                                         geometry,
                                         basis_type,
                                         material,
                                         boundary_points,
                                         internal_points,
                                         positions,
                                         boundary_normal,
                                         surface,
                                         region,
                                         solid_geometry);
        }
    }
}

