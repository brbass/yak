#include "Math_Functions.hh"

#include <cmath>

namespace Math_Functions
{
    using namespace std;

    int factorial(int n)
    {
        return (n == 0) ? 1 : n * factorial(n - 1);
        // int val = 1;
    
        // for (int i = 1; i <= n; ++i)
        // {
        //     val *= i; 
        // }

        // return val;
    }

    double legendre_polynomial(int l,
                               double const &x)
    {
        if (l == 0)
        {
            return 1;
        }
        else if (l == 1)
        {
            return x;
        }
        else
        {
            return ((2 * l - 1) * x * legendre_polynomial(l - 1, x)
                    - (l - 1) * legendre_polynomial(l - 2, x)) / l;
        }
    
        // double plm2 = 0;
        // double plm1 = 1;
        // double pl = x;
    
        // for (int i = 2; i < l + 1; ++i)
        // {
        //     double j = static_cast<double>(i);
        
        //     plm2 = plm1;
        //     plm1 = pl;
        //     pl = ((2 * j - 1) * x * plm1 - (j - 1) * plm2) / j;
        // }
    
        // return pl;
    }

    double legendre_polynomial(int l, int m,
                               double const &x)
    {
        if (m == l)
        {
            return factorial(2 * l) / (pow(2., l) * factorial(l)) * pow(1. - x * x, l / 2.);
        }
        else if (m == l - 1)
        {
            return factorial(2 * l) / (pow(2., l) * factorial(l)) * x * pow(1. - x * x, (l - 1.) / 2.);
        }
        else
        {
            return (2. * (m + 1.) * x / sqrt(1 - x * x) * legendre_polynomial(l, m + 1, x)
                    - legendre_polynomial(l, m + 2)) / ((l - m) * (l + m + 1));
        }
    }

    double spherical_harmonic(int l,
                              int m,
                              double const &mu,
                              double const &phi)
    {
        double val = (m == 0) ? 1 : 0;

        val = sqrt((2 - val) * factorial(l - abs(m)) / factorial(l + abs(m)));

        double t = (m >= 0) ? cos(m * phi) : sin(abs(m) * phi);
    
        return val * legendre_polynomial(l, abs(m), mu) * t;
    }

    // See file sphere_lebedev_rule.cpp for original conversion code
    void xyz_to_spherical(double const &x,
                          double const &y,
                          double const &z,
                          double &mu,
                          double &phi)
    {
        double val = sqrt(x * x + y * y);
        
        if (val < 0)
        {
            phi = acos(x / val);
        }
        else
        {
            phi = acos(x);
        }
        
        if (y < 0)
        {
            phi = -phi;
        }

        mu = z;
    }

    void spherical_to_xyz(double const &mu,
                          double const &phi,
                          double &x,
                          double &y,
                          double &z)
    {
        double val = sqrt(1 - mu * mu);
        
        x = cos(phi) * val;
        y = sin(phi) * val;
        z = mu;
    }

    double spherical_harmonic(int l,
                              int m,
                              double const &x,
                              double const &y,
                              double const &z)
    {
        double mu;
        double phi;

        xyz_to_spherical(x, y, z, mu, phi);
        
        return spherical_harmonic(l, m, mu, phi);
    }
}
