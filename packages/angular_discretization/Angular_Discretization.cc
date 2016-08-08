#include "Angular_Discretization.hh"

#include <cmath>
#include <vector>

#include "Check.hh"
#include "Math_Functions.hh"

using namespace std;

Angular_Discretization::
Angular_Discretization(int dimension, 
                       int number_of_scattering_moments,
                       int number_of_ordinates):
    dimension_(dimension),
    number_of_scattering_moments_(number_of_scattering_moments),
    number_of_ordinates_(number_of_ordinates)
{
    initialize_moment_data();
}

void Angular_Discretization::
initialize_moment_data()
{
    switch(dimension_)
    {
    case 1:
    {
        angular_normalization_ = 2;
        
        int m = 0;
        int sum = 0;
        
        for (int l = 0; l < number_of_scattering_moments_; ++l)
        {
            l_indices_.push_back(l);
            m_indices_.push_back(m);

            sum += 1;
        }
        
        number_of_moments_ = number_of_scattering_moments_;
        
        break;
    }
    case 2:
    {
        angular_normalization_ = 2 * M_PI;
        
        int sum = 0;
        for (int l = 0; l < number_of_scattering_moments_; ++l)
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
        angular_normalization_ = 4 * M_PI;
        
        int sum = 0;
        for (int l = 0; l < number_of_scattering_moments_; ++l)
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
        
        return Math_Functions::legendre_polynomial(mom, mu);
    }
    case 2:
    {
        int l = l_indices_[mom];
        int m = m_indices_[mom];
        
        int o_x = 0 + dimension_ * ord;
        int o_y = 1 + dimension_ * ord;
        double x = ordinates()[o_x];
        double y = ordinates()[o_y];
        double z = sqrt(1 - x * x - y * y);

        return Math_Functions::spherical_harmonic(l, m, x, y, z);
    }
    case 3:
    {
        int l = l_indices_[mom];
        int m = m_indices_[mom];
        
        int o_x = 0 + dimension_ * ord;
        int o_y = 1 + dimension_ * ord;
        int o_z = 2 + dimension_ * ord;
        double x = ordinates()[o_x];
        double y = ordinates()[o_y];
        double z = ordinates()[o_z];
        
        return Math_Functions::spherical_harmonic(l, m, x, y, z);
    }
    default:
        AssertMsg(false, "dimension not found");
        return 0;
    }
}

