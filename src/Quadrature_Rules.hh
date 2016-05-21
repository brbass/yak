#ifndef Quadrature_Rules_hh
#define Quadrature_Rules_hh

#include <vector>

using std::vector;

namespace Quadrature_Rules
{
    /*
      Get vectors of Gauss_Legendre ordinates and weights
    */
    void gauss_legendre(int n, vector<double> &ordinates, vector<double> &weights);
    
    void lebedev(int n, vector<double> &ordinates, vector<double> &weights);
}
#endif


