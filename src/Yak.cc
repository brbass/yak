#include <iostream>
#include <vector>

#include "Gauss_Legendre.hh"
#include "Legendre_Polynomial.hh"

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cerr << "usage: yak [input.xml]" << endl;
        return 1;
    }
    
    // string filename = argv[1];
    
    int val = atoi(argv[1]);

    vector<double> ordinates;
    vector<double> weights;
    
    gauss_legendre_vec(val, ordinates, weights);

    for (int i = 0; i < val; ++i)
    {
        cout << ordinates[i] << "\t" << weights[i] << endl;
    }
    cout << endl;

    double x = 2.1;
    cout << legendre_polynomial(val, x) << endl;
}
