#include "RBF_Mesh.hh"

#include <cmath>
#include <vector>

#include "Check.hh"
#include "Gaussian_RBF.hh"
#include "Inverse_Multiquadric_RBF.hh"
#include "KD_Tree.hh"
#include "Multiquadric_RBF.hh"
#include "Solid_Geometry.hh"
#include "Wendland_RBF.hh"
#include "XML_Functions.hh"

using namespace std;

namespace // anonymous
{
    double get_shape_parameter(int number_of_neighbors,
                               int dimension,
                               double shape_multiplier,
                               vector<double> const &distance)
    {
        // double average_distance = 0;
        
        // for (int i = 0; i < number_of_neighbors; ++i)
        // {
        //     average_distance += sqrt(distance[i]);
        // }
        
        // average_distance /= number_of_neighbors;
        
        // return number_of_neighbors * shape_multiplier / (4 * average_distance);
        
        double average_distance = 0;
        
        int number_to_average = 2 * dimension;
        
        number_to_average = number_to_average > number_of_neighbors ? number_of_neighbors : number_to_average;
        
        for (int i = 0; i < number_to_average; ++i)
        {
            average_distance += sqrt(distance[i]);
        }
        
        average_distance /= number_to_average;
        
        return shape_multiplier / average_distance;
    }

}

RBF_Mesh::
RBF_Mesh(int dimension,
         int number_of_points,
         int number_of_boundary_points,
         int number_of_internal_points,
         int number_of_transition_points,
         int number_of_neighbors,
         int number_of_materials,
         double shape_multiplier,
         Geometry geometry,
         Basis_Type basis_type,
         vector<int> const &material,
         vector<int> const &boundary_points,
         vector<int> const &internal_points,
         vector<int> const &transition_points,
         vector<double> const &positions,
         vector<double> const &boundary_normal,
         vector<double> const &transition_normal,
         vector<int> const &surface,
         vector<int> const &region,
         shared_ptr<Solid_Geometry> const solid_geometry):
    Spatial_Discretization(dimension,
                           geometry),
    number_of_points_(number_of_points),
    number_of_boundary_points_(number_of_boundary_points),
    number_of_internal_points_(number_of_internal_points),
    number_of_transition_points_(number_of_transition_points),
    number_of_neighbors_(number_of_neighbors),
    number_of_materials_(number_of_materials),
    shape_multiplier_(shape_multiplier),
    material_(material),
    boundary_points_(boundary_points),
    internal_points_(internal_points),
    transition_points_(transition_points),
    surface_(surface),
    region_(region),
    point_positions_(positions),
    boundary_normal_(boundary_normal),
    transition_normal_(transition_normal),
    boundary_nodes_(number_of_boundary_points, true),
    solid_geometry_(solid_geometry)
{
    Check(0 < number_of_neighbors <= number_of_points);

    
    kd_tree_ = make_shared<KD_Tree>(dimension,
                                    number_of_points,
                                    positions);
    
    neighbors_.resize(number_of_points_);
    shape_parameter_.resize(number_of_points_ * dimension_);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<int> neighbors(number_of_neighbors_);
        vector<double> distances(number_of_neighbors_);
        
        kd_tree_->find_neighbors(i,
                                 number_of_neighbors_,
                                 neighbors,
                                 distances);
        
        neighbors_[i] = neighbors;
        
        double shape_parameter = get_shape_parameter(number_of_neighbors_,
                                                     dimension_,
                                                     shape_multiplier_,
                                                     distances);
        
        for (int d = 0; d < dimension_; ++d)
        {
            int k = d + dimension_ * i;
            
            shape_parameter_[k] = shape_parameter;
        }
    }
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        vector<double> position(dimension_);
        vector<double> shape(dimension_);
        for (int d = 0; d < dimension_; ++d)
        {
            position[d] = positions[d + dimension_ * i];
            shape[d] = shape_parameter_[d + dimension_ * i];
        }
        
        shared_ptr<RBF> rbf;
        
        switch(basis_type)
        {
        case Basis_Type::GAUSSIAN:
        {
            rbf = make_shared<Gaussian_RBF>(dimension_,
                                            position,
                                            shape);
            break;
        }
        case Basis_Type::MULTIQUADRIC:
        {
            rbf = make_shared<Multiquadric_RBF>(dimension,
                                                position,
                                                shape);
            break;
        }
        case Basis_Type::INVERSE_MULTIQUADRIC:
        {
            rbf = make_shared<Inverse_Multiquadric_RBF>(dimension,
                                                        position,
                                                        shape);
            break;
        }
        case Basis_Type::WENDLAND30:
        {
            int order = 0;
            
            rbf = make_shared<Wendland_RBF>(dimension,
                                            order,
                                            position,
                                            shape);
            break;
        }
        case Basis_Type::WENDLAND31:
        {
            int order = 1;
            
            rbf = make_shared<Wendland_RBF>(dimension,
                                            order,
                                            position,
                                            shape);
            break;
        }
        case Basis_Type::WENDLAND32:
        {
            int order = 2;
            
            rbf = make_shared<Wendland_RBF>(dimension,
                                            order,
                                            position,
                                            shape);
            break;
        }
        case Basis_Type::WENDLAND33:
        {
            int order = 3;
            
            rbf = make_shared<Wendland_RBF>(dimension,
                                            order,
                                            position,
                                            shape);
            break;
        }
        default:
            AssertMsg(false, "no such type of basis function: ");
            break;
        }

        basis_functions_.push_back(rbf);
    }
    
    check_class_invariants();

    get_cell_type(cell_type_);
}


RBF_Mesh::
RBF_Mesh(int dimension,
         int number_of_points,
         int number_of_boundary_points,
         int number_of_internal_points,
         int number_of_neighbors,
         int number_of_materials,
         double shape_multiplier,
         Geometry geometry,
         Basis_Type basis_type,
         vector<int> const &material,
         vector<int> const &boundary_points,
         vector<int> const &internal_points,
         vector<double> const &positions,
         vector<double> const &boundary_normal,
         vector<int> const &surface,
         vector<int> const &region,
         shared_ptr<Solid_Geometry> const solid_geometry):
    RBF_Mesh(dimension,
             number_of_points,
             number_of_boundary_points,
             number_of_internal_points,
             0, // number_of_transition_points
             number_of_neighbors,
             number_of_materials,
             shape_multiplier,
             geometry,
             basis_type,
             material,
             boundary_points,
             internal_points,
             vector<int>(), // transition_points
             positions,
             boundary_normal,
             vector<double>(), // transition_normal
             surface,
             region,
             solid_geometry)
{
}

void RBF_Mesh::
check_class_invariants() const
{
    Assert(number_of_boundary_points_ + number_of_internal_points_ + number_of_transition_points_ == number_of_points_);
    Assert(boundary_nodes_.size() == number_of_boundary_points_);
    Assert(boundary_points_.size() == number_of_boundary_points_);
    Assert(internal_points_.size() == number_of_internal_points_);
    Assert(transition_points_.size() == number_of_transition_points_);
    Assert(material_.size() == number_of_points_);
    Assert(boundary_normal_.size() == number_of_boundary_points_ * dimension_);
    Assert(transition_normal_.size() == number_of_transition_points_ * dimension_);
    Assert(point_positions_.size() == number_of_points_ * dimension_);
    Assert(shape_parameter_.size() == number_of_points_ * dimension_);
    Assert(basis_functions_.size() == number_of_points_);
}

void RBF_Mesh::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node rbf = output_node.append_child("rbf_mesh");
    
    XML_Functions::append_child(rbf, dimension_, "dimension");
    XML_Functions::append_child(rbf, number_of_points_, "number_of_points");
    XML_Functions::append_child(rbf, number_of_boundary_points_, "number_of_boundary_points");
    XML_Functions::append_child(rbf, number_of_internal_points_, "number_of_internal_points");
    XML_Functions::append_child(rbf, number_of_transition_points_, "number_of_transition_points");
    XML_Functions::append_child(rbf, number_of_materials_, "number_of_materials");
    XML_Functions::append_child(rbf, boundary_points_, "boundary_points", "point");
    XML_Functions::append_child(rbf, internal_points_, "internal_points", "point");
    XML_Functions::append_child(rbf, transition_points_, "transition_points", "point");
    XML_Functions::append_child(rbf, material_, "material", "point");
    XML_Functions::append_child(rbf, boundary_normal_, "boundary_normal", "point-dimension");
    XML_Functions::append_child(rbf, transition_normal_, "transition_normal", "point-dimension");
    XML_Functions::append_child(rbf, point_positions_, "point_positions", "dimension-point");
    XML_Functions::append_child(rbf, shape_parameter_, "shape_parameter", "dimension-point");
}

void RBF_Mesh::
set_shape_multiplier(double shape_multiplier)
{
    double old_shape_multiplier = shape_multiplier_;
    shape_multiplier_ = shape_multiplier;

    double mult = shape_multiplier / old_shape_multiplier;
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        for (int d = 0; d < dimension_; ++d)
        {
            int k_sm = d + dimension_ * i;
            
            shape_parameter_[k_sm] = shape_parameter_[k_sm] * mult;
        }
    }

    for (int i = 0; i < number_of_points_; ++i)
    {
        basis_functions_[i]->multiply_shape_parameter(mult);
    }
}
