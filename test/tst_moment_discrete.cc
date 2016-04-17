#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Discrete_To_Moment.hh"
#include "Energy_Discretization.hh"
#include "Finite_Element_Mesh.hh"
#include "Gauss_Legendre_Ordinates.hh"
#include "Moment_To_Discrete.hh"
#include "Random_Number_Generator.hh"

using namespace std;

void print(vector<double> &x0, vector<double> &x)
{
    for (unsigned i = 0; i < x.size(); ++i)
    {
        cout << setw(16) << x0[i];
        cout << setw(16) << x[i];
        cout << endl;
    }
    
    cout << endl;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cerr << "usage: tst_moment_discrete [num_moments num_ordinates]" << endl;
        return 1;
    }
    
    // string filename = argv[1];
    
    int dimension = 1;
    int number_of_cells = 10;
    int number_of_nodes = 2;
    int number_of_groups = 2;
    int number_of_moments = atoi(argv[1]);
    int number_of_ordinates = atoi(argv[2]);
    double length = 1;

    Random_Number_Generator rng(0, 1);
    vector<double> x0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_moments);
    vector<double> y0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_ordinates);
    
    shared_ptr<Spatial_Discretization> spatial_discretization 
        = make_shared<Finite_Element_Mesh>(dimension,
                                           number_of_cells,
                                           number_of_nodes,
                                           length,
                                           "DFEM");
    shared_ptr<Angular_Discretization> angular_discretization 
        = make_shared<Gauss_Legendre_Ordinates>(dimension,
                                                number_of_moments,
                                                number_of_ordinates);
    shared_ptr<Energy_Discretization> energy_discretization
        = make_shared<Energy_Discretization>(number_of_groups);
    
    shared_ptr<Vector_Operator> upM 
        = make_shared<Moment_To_Discrete>(spatial_discretization,
                                          angular_discretization,
                                          energy_discretization);
    
    shared_ptr<Vector_Operator> upD 
        = make_shared<Discrete_To_Moment>(spatial_discretization,
                                          angular_discretization,
                                          energy_discretization);
    
    vector<double> x(x0);
    
    (*upD)(x);
    (*upM)(x);
    print(x0, x);
}
