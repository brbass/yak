#include "RBF_Sweep_1D.hh"

#include "Check.hh"
#include "Dense_Solve.hh"

using namespace std;

RBF_Sweep_1D::
RBF_Sweep_1D(shared_ptr<Spatial_Discretization> spatial_discretization,
             shared_ptr<Angular_Discretization> angular_discretization,
             shared_ptr<Energy_Discretization> energy_discretization,
             shared_ptr<Nuclear_Data> nuclear_data,
             shared_ptr<Source_Data> source_data):
    Vector_Operator(get_size(spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             source_data),
                    get_size(spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             source_data)),
    Ordinate_Sweep_Operator(spatial_discretization,
                            angular_discretization,
                            energy_discretization,
                            nuclear_data,
                            source_data)
{
    rbf_mesh_ = dynamic_pointer_cast<RBF_Mesh>(spatial_discretization);
    Assert(rbf_mesh_);
}

void RBF_Sweep_1D::
apply(vector<double> &x)
{
    switch(spatial_discretization_->geometry())
    {
    case Spatial_Discretization::Geometry::SLAB:
        sweep_slab(x);
        break;
    default:
        AssertMsg(false, "Sweep type not implemented");
        break;
    }
}

void RBF_Sweep_1D::
sweep_slab(vector<double> &x)
{
    int number_of_points = spatial_discretization_->number_of_cells();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_boundary_points = spatial_discretization_->number_of_boundary_cells();
    int number_of_augments = source_data_->number_of_augments();
    int psi_size = row_size() - number_of_augments;
    vector<double> const ordinates = angular_discretization_->ordinates();
    vector<double> const sigma_t = nuclear_data_->sigma_t();
    vector<double> const boundary_source = source_data_->boundary_source();
    vector<double> const alpha = source_data_->alpha();
    
    Dense_Solve matrix_solver(number_of_points);
    
    int d = 0; // dimension
    vector<double> a_data(number_of_points * number_of_points, 0.0); // matrix
    vector<double> b_data(number_of_points, 0.0); // rhs
    vector<double> x_data(number_of_points, 1.0); //lhs
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            // a_data.assign(number_of_points * number_of_points, 0.0);
            
            if (ordinates[o]>0)
            {
                for (int i = 0; i < number_of_points; ++i) // basis point
                {
                    shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(i);
                    
                    // internal points
                    for (int j = 1; j < number_of_points; ++j) // equation point
                    {
                        shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(j);
                        vector<double> const equation_position = equation_rbf->position();
                        
                        int k_a = i + number_of_points * j;
                        int k_sig = g + number_of_groups * j;

                        a_data[k_a] = ordinates[o] * basis_rbf->dbasis(d, equation_position)
                            + sigma_t[k_sig] * basis_rbf->basis(equation_position);
                    }
                    
                    // boundary point
                    int j = 0;
                    int k_a = i;
                    
                    shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(j);
                    vector<double> const equation_position = equation_rbf->position();
                    
                    a_data[k_a] = basis_rbf->basis(equation_position);
                }
            }
            else
            {
                for (int i = 0; i < number_of_points; ++i) // basis point
                {
                    shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(i);
                    
                    // internal points
                    for (int j = 0; j < number_of_points - 1; ++j) // equation point
                    {
                        shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(j);
                        vector<double> const equation_position = equation_rbf->position();

                        int k_a = i + number_of_points * j;
                        int k_sig = g + number_of_groups * j;
                        
                        a_data[k_a] = ordinates[o] * basis_rbf->dbasis(d, equation_position)
                            + sigma_t[k_sig] * basis_rbf->basis(equation_position);
                    }
                    
                    // boundary point
                    int j = number_of_points - 1;
                    int k_a = i + number_of_points * j;
                    
                    shared_ptr<RBF> equation_rbf = rbf_mesh_->basis_function(j);
                    vector<double> const equation_position = equation_rbf->position();
                    
                    a_data[k_a] = basis_rbf->basis(equation_position);
                }
            }
            
            // initialize rhs

            if (ordinates[o]>0)
            {
                // internal points
                for (int j=1; j<number_of_points; ++j) // equation point
                {
                    int k_x = g + number_of_groups * (o + number_of_ordinates * j);
                    
                    b_data[j] = x[k_x];
                }
                
                // boundary point
                int j = 0;
                int o1 = number_of_ordinates - o - 1;
                int b = 0;
                int k_ref = psi_size + g + number_of_groups * (o1 + number_of_ordinates * b);
                b_data[j] = alpha[b] * x[k_ref];
                
                if (include_boundary_source_)
                {
                    int k_bs = g + number_of_groups * (o + number_of_ordinates * b);
                    
                    b_data[j] += boundary_source[k_bs];
                }
            }
            else
            {
                // internal points
                for (int j=0; j<number_of_points-1; ++j)
                {
                    int k_x = g + number_of_groups * (o + number_of_ordinates * j);
                    
                    b_data[j] = x[k_x];
                }
                
                // boundary point
                int j = number_of_points - 1;
                int o1 = number_of_ordinates - o - 1;
                int b = 1;
                int k_ref = psi_size + g + number_of_groups * (o1 + number_of_ordinates * b);
                
                b_data[j] = alpha[b] * x[k_ref];
                
                if (include_boundary_source_)
                {
                    int k_bs = g + number_of_groups * (o + number_of_ordinates * b);
                    
                    b_data[j] += boundary_source[k_bs];
                }
            }
            
            // initialize lhs

            x_data.assign(number_of_points, 1.0);
            
            // perform matrix solve
            
            matrix_solver.solve(a_data, b_data, x_data);

            for (int i = 0; i < number_of_points; ++i)
            {
                shared_ptr<RBF> location_rbf = rbf_mesh_->basis_function(i);
                vector<double> const position = location_rbf->position();
                
                double sum = 0;
                
                for (int j = 0; j < number_of_points; ++j)
                {
                    shared_ptr<RBF> basis_rbf = rbf_mesh_->basis_function(j);
                    
                    sum += x_data[j] * basis_rbf->basis(position);
                }
                
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                
                x[k_x] = sum;
            }
        }
    }

    // update augments

    vector<int> boundary_points = spatial_discretization_->boundary_cells();
    for (int b = 0; b < number_of_boundary_points; ++b)
    {
        int i = boundary_points[b];
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_b = psi_size + g + number_of_groups * (o + number_of_ordinates * b);
                int k_psi = g + number_of_groups * (o + number_of_ordinates * i);
                    
                x[k_b] = x[k_psi];
            }
        }
    }
}
