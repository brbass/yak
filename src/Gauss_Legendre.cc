#include "Gauss_Legendre.hh"

#include <iostream>
#include <memory>
#include <vector>

#include "gauss_legendre.hpp"

#include "Check.hh"

void gauss_legendre_vec(int n, std::vector<double> &ordinates, std::vector<double> &weights)
{
    if (n % 2 == 0)
    {
        double *x = new double[n/2];
        double *w = new double[n/2];
        double eps = 1e-15;
        
        gauss_legendre_tbl(n, x, w, eps);
        ordinates.resize(n, 0);
        weights.resize(n, 0);
    
        for (unsigned o = 0; o < n / 2; ++o)
        {
            ordinates[o + n/2] = x[o];
            weights[o + n/2] = w[o];
            ordinates[o] = -x[n/2 - o - 1];
            weights[o] = w[n/2 - o - 1];
        }
    
        delete[] x;
        delete[] w;
    }
    else
    {
        AssertMsg(false, "Only even quadratures supported");
    }
}

