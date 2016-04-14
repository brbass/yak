#include "Legendre_Polynomial.hh"

double legendre_polynomial(int n, double &x)
{
    if (n == 0)
    {
        return 1;
    }
    else if (n == 1)
    {
        return x;
    }
    
    double pnm2 = 0;
    double pnm1 = 1;
    double pn = x;
    
    for (int i = 2; i < n + 1; ++i)
    {
        double j = static_cast<double>(i);
        
        pnm2 = pnm1;
        pnm1 = pn;
        pn = ((2 * j - 1) * x * pnm1 - (j - 1) * pnm2) / j;
    }
    
    return pn;
}
