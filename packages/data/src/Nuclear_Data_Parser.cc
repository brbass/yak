#include "Nuclear_Data_Parser.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Nuclear_Data_Parser::
Nuclear_Data_Parser(pugi::xml_node &input_file,
                    shared_ptr<Spatial_Discretization> spatial,
                    shared_ptr<Angular_Discretization> angular,
                    shared_ptr<Energy_Discretization> energy):
    Parser(input_file),
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    pugi::xml_node nuclear_node = input_file.child("nuclear_data");
    pugi::xml_node materials = nuclear_node.child("materials");

    int number_of_materials = XML_Functions::child_value<int>(materials, "number_of_materials");
    Assert(number_of_materials = spatial_->number_of_materials());
    
    int number_of_cells = spatial_->number_of_cells();
    int number_of_transition_cells = spatial_->number_of_transition_points();
    int number_of_moments = angular_->number_of_scattering_moments();
    int number_of_groups = energy_->number_of_groups();
    
    // parse the data for each material

    vector<double> sigma_t_m(number_of_materials * number_of_groups);
    vector<double> sigma_s_m(number_of_materials * number_of_groups * number_of_groups * number_of_moments);
    vector<double> nu_m(number_of_materials * number_of_groups);
    vector<double> sigma_f_m(number_of_materials * number_of_groups);
    vector<double> chi_m(number_of_materials * number_of_groups);
    
    for (pugi::xml_node material = materials.child("material"); material; material = material.next_sibling("material"))
    {
        int a = XML_Functions::child_value<int>(material, "material_number");
        
        vector<double> sigma_t_a = XML_Functions::child_vector<double>(material, "sigma_t", number_of_groups);
        vector<double> sigma_s_a = XML_Functions::child_vector<double>(material, "sigma_s", number_of_groups * number_of_groups * number_of_moments);
        vector<double> nu_a = XML_Functions::child_vector<double>(material, "nu", number_of_groups);
        vector<double> sigma_f_a = XML_Functions::child_vector<double>(material, "sigma_f", number_of_groups);
        vector<double> chi_a = XML_Functions::child_vector<double>(material, "chi", number_of_groups);
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * a;
            
            sigma_t_m[k] = sigma_t_a[g];
            nu_m[k] = nu_a[g];
            sigma_f_m[k] = sigma_f_a[g];
            chi_m[k] = chi_a[g];
        }
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int gf = 0; gf < number_of_groups; ++gf)
                {
                    int k_a = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * a));
                    int k_l = gf + number_of_groups * (gt + number_of_groups * m);
                    
                    sigma_s_m[k_a] = sigma_s_a[k_l];
                }
            }
        }
    } // materials
    
    // assign the data for each cell from the material data
    
    vector<int> const material = spatial_->material();

    int cells_plus_transition = number_of_cells + number_of_transition_cells;
    vector<double> sigma_t(cells_plus_transition * number_of_groups);
    vector<double> sigma_s(cells_plus_transition * number_of_groups * number_of_groups * number_of_moments);
    vector<double> nu(cells_plus_transition * number_of_groups);
    vector<double> sigma_f(cells_plus_transition * number_of_groups);
    vector<double> chi(cells_plus_transition * number_of_groups);

    int number_of_internal_cells = spatial_->number_of_internal_points();
    int number_of_boundary_cells = spatial_->number_of_boundary_points();
    vector<int> const internal_cells = spatial_->internal_cells();
    vector<int> const boundary_cells = spatial_->boundary_cells();
    
    for (int c = 0; c < number_of_internal_cells; ++c)
    {
        int i = internal_cells[c];
        int a = material[i];
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k_a = g + number_of_groups * a;
            int k_i = g + number_of_groups * i;
            
            sigma_t[k_i] = sigma_t_m[k_a];
            nu[k_i] = nu_m[k_a];
            sigma_f[k_i] = sigma_f_m[k_a];
            chi[k_i] = chi_m[k_a];
        }
        
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int gf = 0; gf < number_of_groups; ++gf)
                {
                    int k_a = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * a));
                    int k_i = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * i));
                    
                    sigma_s[k_i] = sigma_s_m[k_a];
                }
            }
        }
    } // internal cells
    
    for (int c = 0; c < number_of_boundary_cells; ++c)
    {
        int i = boundary_cells[c];
        int a = material[i];
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k_a = g + number_of_groups * a;
            int k_i = g + number_of_groups * i;
            
            sigma_t[k_i] = sigma_t_m[k_a];
            nu[k_i] = nu_m[k_a];
            sigma_f[k_i] = sigma_f_m[k_a];
            chi[k_i] = chi_m[k_a];
        }
        
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int gf = 0; gf < number_of_groups; ++gf)
                {
                    int k_a = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * a));
                    int k_i = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * i));
                    
                    sigma_s[k_i] = sigma_s_m[k_a];
                }
            }
        }
    } // boundary cells
    
    vector<int> transition_cells = spatial_->transition_cells();
    
    for (int c = 0; c < number_of_transition_cells; ++c)
    {
        int i = transition_cells[c];
        int a = material[i];
        
        int a2 = a % number_of_materials;
        int a1 = a - number_of_materials * a2;
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k_a1 = g + number_of_groups * a2;
            int k_i1 = g + number_of_groups * i;
            int k_a2 = g + number_of_groups * a2;
            int k_i2 = g + number_of_groups * (number_of_cells + c);
            
            sigma_t[k_i1] = sigma_t_m[k_a1];
            nu[k_i1] = nu_m[k_a1];
            sigma_f[k_i1] = sigma_f_m[k_a1];
            chi[k_i1] = chi_m[k_a1];

            sigma_t[k_i2] = sigma_t_m[k_a2];
            nu[k_i2] = nu_m[k_a2];
            sigma_f[k_i2] = sigma_f_m[k_a2];
            chi[k_i2] = chi_m[k_a2];
        }
        
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int gf = 0; gf < number_of_groups; ++gf)
                {
                    int k_a1 = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * a1));
                    int k_i1 = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * i));
                    int k_a2 = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * a2));
                    int k_i2 = gf + number_of_groups * (gt + number_of_groups * (m + number_of_moments * (number_of_cells + c)));
                    
                    sigma_s[k_i1] = sigma_s_m[k_a1];
                    sigma_s[k_i2] = sigma_s_m[k_a2];
                }
            }
        }
    }
    
    nuclear_ = make_shared<Nuclear_Data>(spatial_,
                                         angular_,
                                         energy_,
                                         sigma_t,
                                         sigma_s,
                                         nu,
                                         sigma_f,
                                         chi);
} // constructor
