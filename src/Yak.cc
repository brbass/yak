#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Discrete_To_Moment.hh"
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
        cerr << "usage: yak [input.xml]" << endl;
        return 1;
    }
    
    // string filename = argv[1];
    
    int number_of_cells = 10;
    int number_of_nodes = 2;
    int number_of_groups = 2;
    int number_of_moments = atoi(argv[1]);
    int number_of_ordinates = atoi(argv[2]);

    Random_Number_Generator rng(0, 1);
    vector<double> x0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_moments);
    vector<double> y0 = rng.random_double_vector(number_of_cells * number_of_nodes * number_of_groups * number_of_ordinates);
    
    shared_ptr<Vector_Operator> upM 
        = make_shared<Moment_To_Discrete>(number_of_cells,
                                          number_of_nodes,
                                          number_of_groups,
                                          number_of_moments,
                                          number_of_ordinates);
    
    shared_ptr<Vector_Operator> upD
        = make_shared<Discrete_To_Moment>(number_of_cells,
                                          number_of_nodes,
                                          number_of_groups,
                                          number_of_moments,
                                          number_of_ordinates);
    
    vector<double> x(x0);
    
    (*upM)(x);
    (*upD)(x);
    print(x0, x);
}
