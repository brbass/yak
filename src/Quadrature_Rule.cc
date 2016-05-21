#include "Quadrature_Rule.hh"

#include <iostream>
#include <memory>
#include <vector>

#include "quadrule.hpp"
#include "sphere_lebedev_rule.hpp"

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
            legendre_set(n, &ordinates[0], &weights[0]);
        }
        else
        {
            legendre_dr_compute(n, &ordinates[0], &weights[0]);
        }
    }

    void lebedev(int n,
                 int dimension,
                 int &number_of_ordinates,
                 vector<double> &ordinates,
                 vector<double> &weights)
    {
        Assert(available_table(n) == 1);

        number_of_ordinates = order_table(n);

        ordinates.resize(number_of_ordinates * dimension);
        weights.resize(number_of_ordinates);
        
        vector<double> x(number_of_ordinates);
        vector<double> y(number_of_ordinates);
        vector<double> z(number_of_ordinates);

        ld_by_order(number_of_ordinates, &x[0], &y[0], &z[0], &weights[0]);

        if (dimension == 1)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                ordinates[o] = x[o];
            }
        }
        else if (dimension == 2)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                ordinates[0 + dimension * o] = x[o];
                ordinates[1 + dimension * o] = x[o];
            }
        }
        else if (dimension == 3)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                ordinates[0 + dimension * o] = x[o];
                ordinates[1 + dimension * o] = x[o];
                ordinates[2 + dimension * o] = x[o];
            }
        }
    }
}
