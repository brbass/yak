#ifndef Legendre_Polynomial_hh
#define Legendre_Polynomial_hh

// Factorial
int factorial(int n);

// Returns Legendre polynomial
double legendre_polynomial(int l, double &x);

// Returns associated real Legendre polynomial function
double legendre_polynomial(int l, int m, double &x);

// Returns real spherical harmonic function
double spherical_harmonic(int l, int m, double &mu, double &phi);

#endif
