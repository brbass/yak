#include "Source_Data_Parser.hh"

#include "RBF_Mesh.hh"
#include "Solid_Geometry.hh"

using namespace std;

Source_Data_Parser::
Source_Data_Parser(pugi::xml_node &input_file,
                    shared_ptr<Spatial_Discretization> spatial,
                    shared_ptr<Angular_Discretization> angular,
                    shared_ptr<Energy_Discretization> energy):
    Parser(input_file),
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    pugi::xml_node source_node = input_file.child("source_data");
    pugi::xml_node nuclear_node = input_file.child("nuclear_data");
    pugi::xml_node materials = nuclear_node.child("materials");

    int number_of_materials = XML_Functions::child_value<int>(materials, "number_of_materials");
    int dimension = spatial_->dimension();
    int number_of_cells = spatial_->number_of_cells();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_boundary_cells = spatial_->number_of_boundary_cells();
    int number_of_moments = angular_->number_of_moments();
    int number_of_ordinates = angular_->number_of_ordinates();
    int number_of_groups = energy_->number_of_groups();
    vector<bool> const boundary_nodes = spatial_->boundary_nodes();
    
    string internal_type_str = XML_Functions::child_value<string>(source_node, "internal_source_type");
    string boundary_type_str = XML_Functions::child_value<string>(source_node, "boundary_source_type");
    
    Source_Data::Source_Type internal_type;
    Source_Data::Source_Type boundary_type;
    
    vector<double> internal_source;

    if (internal_type_str == "isotropic_full")
    {
        internal_type = Source_Data::Source_Type::FULL;
        
        double angular_normalization = angular_->angular_normalization();
        
        vector<double> internal_m(number_of_materials * number_of_groups);
        
        for (pugi::xml_node material = materials.child("material"); material; material = material.next_sibling("material"))
        {
            int a = XML_Functions::child_value<int>(material, "material_number");
            
            vector<double> internal_a = XML_Functions::child_vector<double>(material, "internal_source", number_of_groups);
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * a;
                
                internal_m[k] = internal_a[g] / angular_normalization;
            }
        }

        internal_source.resize(number_of_cells * number_of_ordinates * number_of_groups * number_of_nodes);
        
        vector<int> const material = spatial_->material();

        for (int i = 0; i < number_of_cells; ++i)
        {
            int a = material[i];
            
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_m = g + number_of_groups * a;
                    
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k_i = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                        
                        internal_source[k_i] = internal_m[k_m];
                    }
                }
            }
        }
    }
    else if (internal_type_str == "isotropic")
    {
        internal_type = Source_Data::Source_Type::MOMENT;
        
        // double angular_normalization = angular_->angular_normalization();
        
        vector<double> internal_m(number_of_materials * number_of_groups);
        
        for (pugi::xml_node material = materials.child("material"); material; material = material.next_sibling("material"))
        {
            int a = XML_Functions::child_value<int>(material, "material_number");
            
            vector<double> internal_a = XML_Functions::child_vector<double>(material, "internal_source", number_of_groups);
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * a;
                
                internal_m[k] = internal_a[g];
            }
        }

        internal_source.assign(number_of_cells * number_of_moments * number_of_groups * number_of_nodes, 0);
        
        vector<int> const material = spatial_->material();

        for (int i = 0; i < number_of_cells; ++i)
        {
            int a = material[i];
            
            {
                int m = 0;
                
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_m = g + number_of_groups * a;
                    
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k_i = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                        
                        internal_source[k_i] = internal_m[k_m];
                    }
                }
            }
        }
    }
    else 
    {
        AssertMsg(false, "source type \"" + internal_type_str + "\" not found");
    }
    
    vector<double> boundary_source;
    vector<double> alpha;

    if (boundary_type_str == "cellwise_isotropic")
    {
        boundary_type = Source_Data::Source_Type::FULL;
        
        // double angular_normalization = angular_->angular_normalization();
        
        vector<double> boundary_m = XML_Functions::child_vector<double>(source_node, "boundary_source", number_of_boundary_cells * number_of_groups);
        
        boundary_source.assign(number_of_boundary_cells * number_of_nodes * number_of_ordinates * number_of_groups, 0);
        
        for (int b = 0; b < number_of_boundary_cells; ++b)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                int k_bn = n + number_of_nodes * b;
                
                if (boundary_nodes[k_bn])
                {
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        for (int g = 0; g < number_of_groups; ++g)
                        {
                            int k = g + number_of_groups * b;
                            int k_o = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * b));
                    
                            boundary_source[k_o] = boundary_m[k];
                        }
                    }
                }
            }
        }
    }
    else if (boundary_type_str == "surface_isotropic")
    {
        boundary_type = Source_Data::Source_Type::FULL;

        vector<int> const boundary_cells = spatial->boundary_cells();
        shared_ptr<RBF_Mesh> const rbf_mesh = dynamic_pointer_cast<RBF_Mesh>(spatial);
        shared_ptr<Solid_Geometry> const solid_geometry = rbf_mesh->solid_geometry();
        int number_of_surfaces = solid_geometry->number_of_surfaces();
        vector<int> const surfaces = rbf_mesh->surface();

        vector<double> boundary_s = XML_Functions::child_vector<double>(source_node, "boundary_source", number_of_surfaces * number_of_groups);

        boundary_source.assign(number_of_boundary_cells * number_of_nodes * number_of_ordinates * number_of_groups, 0);
        alpha.assign(number_of_boundary_cells, 0);
        
        for (int b = 0; b < number_of_boundary_cells; ++b)
        {
            int i = boundary_cells[b];
            shared_ptr<RBF> const basis = rbf_mesh->basis_function(i);
            int surface = surfaces[b];
            
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k = g + number_of_groups * surface;
                    int k_o = g + number_of_groups * (o + number_of_ordinates * b);
                    
                    boundary_source[k_o] = boundary_s[k];
                }
            }

            Surface::Surface_Type surface_type = solid_geometry->surface(surface)->surface_type();

            if (surface_type == Surface::Surface_Type::REFLECTIVE_BOUNDARY)
            {
                alpha[b] = 1.0;
            }
            else
            {
                alpha[b] = 0.0;
            }
        }
    }
    else 
    {
        AssertMsg(false, "source type \"" + boundary_type_str + "\" not found");
    }
    
    if (dimension == 1)
    {
        alpha = XML_Functions::child_vector<double>(source_node, "alpha", number_of_boundary_cells);
    }
    
    source_ = make_shared<Source_Data>(internal_type,
                                       boundary_type,
                                       spatial_,
                                       angular_,
                                       energy_,
                                       internal_source,
                                       boundary_source,
                                       alpha);
}
