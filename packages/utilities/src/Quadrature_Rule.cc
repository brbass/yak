#include "Quadrature_Rule.hh"

#include <iostream>
#include <memory>
#include <vector>

#include "quadrule.hh"
#include "sphere_lebedev_rule.hh"

#include "Check.hh"

namespace Quadrature_Rule
{
    void gauss_legendre(int n,
                        vector<double> &ordinates,
                        vector<double> &weights)
    {
        if (n < 1)
        {
            AssertMsg(false, "quadrature must have n >= 1");
        }

        ordinates.resize(n);
        weights.resize(n);
    
        if (n <= 33
            || (63 <= n && n <= 65)
            || (127 <= n && n <= 129)
            || (255 <= n && 257 <= n))
        {
            quadrule::legendre_set(n, &ordinates[0], &weights[0]);
        }
        else
        {
            quadrule::legendre_dr_compute(n, &ordinates[0], &weights[0]);
        }
    }

    int lebedev(int n,
                int dimension,
                vector<double> &ordinates,
                vector<double> &weights)
    {
        if (n < 1)
        {
            AssertMsg(false, "quadrature must have n >= 1");
        }
        if (dimension != 2 && dimension != 3)
        {
            AssertMsg(false, "Lebedev quadrature must have d = 2 or 3");
        }
        
        Assert(sphere_lebedev_rule::available_table(n) == 1);

        int num = sphere_lebedev_rule::order_table(n);
        int number_of_ordinates = num;
        
        ordinates.resize(0);
        weights.resize(0);
        
        vector<double> x(number_of_ordinates);
        vector<double> y(number_of_ordinates);
        vector<double> z(number_of_ordinates);
        vector<double> w(number_of_ordinates);
        
        sphere_lebedev_rule::ld_by_order(number_of_ordinates, &x[0], &y[0], &z[0], &w[0]);

        if (dimension == 2)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                if (z[o] >= 0)
                {
                    ordinates.push_back(x[o]);
                    ordinates.push_back(y[o]);
                    weights.push_back(w[o]);
                }
            }
            number_of_ordinates = weights.size();
        }
        else // (dimension == 3)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                ordinates.push_back(x[o]);
                ordinates.push_back(y[o]);
                ordinates.push_back(z[o]);
                weights.push_back(w[o]);
            }
        }

        return number_of_ordinates;
    }
}
