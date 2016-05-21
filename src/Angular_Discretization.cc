#include "Angular_Discretization.hh"

#include <cmath>
#include <vector>

#include "Check.hh"
#include "Legendre_Polynomial.hh"

using namespace std;

Angular_Discretization::
Angular_Discretization(int dimension, 
                       int number_of_scattering_moments,
                       int number_of_ordinates):
    dimension_(dimension),
    number_of_scattering_moments_(number_of_scattering_moments),
    number_of_ordinates_(number_of_ordinates)
{
    switch(dimension_)
    {
    case 1:
    {
        angular_normalization_ = 2;

        int m = 0;
        int sum = 0;
        
        for (int l = 0; l < number_of_scattering_moments; ++l)
        {
            l_indices_.push_back(l);
            m_indices_.push_back(m);

            sum += 1;
        }

        number_of_moments_ = sum;
        
        break;
    }
    case 2:
    {
        int sum = 0;
        for (int l = 0; l < number_of_scattering_moments; ++l)
        {
            for (int m = 0; m <= l; ++m)
            {
                l_indices_.push_back(l);
                m_indices_.push_back(m);

                sum += 1;
            }
        }
        
        number_of_moments_ = sum;
        
        break;
    }
    case 3:
    {
        int sum = 0;
        for (int l = 0; l < number_of_scattering_moments; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                l_indices_.push_back(l);
                m_indices_.push_back(m);
                
                sum += 1;
            }
        }
        
        number_of_moments_ = sum;
        
        break;
    }
    default:
        AssertMsg(false, "Dimension not found");
        break;
    }
}

double Angular_Discretization::
moment(int mom,
       int ord)
{
    switch(dimension_)
    {
    case 1:
    {
        double mu = ordinates()[ord];
        
        return legendre_polynomial(mom, mu);
    }
    case 2:
    {
        int l = l_indices[mom];
        int m = m_indices[mom];
        
        int o_x = 0 + dimension_ * ord;
        int o_y = 1 + dimension_ * ord;
        double mu = ordinates()[o_x];
        double phi = acos(ordinates()[o_y] / sqrt(1 - mu * mu));
        
        return spherical_harmonic(l, m, mu, phi);
    }
    case 3:
    {
        int l = l_indices[mom];
        int m = m_indices[mom];
        
        int o_x = 0 + dimension_ * ord;
        int o_y = 1 + dimension_ * ord;
        double mu = ordinates()[o_x];
        double phi = acos(ordinates()[o_y] / sqrt(1 - mu * mu));
        
        return spherical_harmonic(l, m, mu, phi);
    }
    default:
        AssertMsg(false, "dimension not found");
        return 0;
    }
}

