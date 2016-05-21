#ifndef Math_Functions_hh
#define Math_Functions_hh

namespace Math_Functions
{
    // Factorial
    int factorial(int n);
    
    // Returns Legendre polynomial
    double legendre_polynomial(int l, double &x);

    // Returns associated real Legendre polynomial function
    double legendre_polynomial(int l, int m, double &x);

    // Returns real spherical harmonic function
    double spherical_harmonic(int l, int m, double &mu, double &phi);
}
#endif
