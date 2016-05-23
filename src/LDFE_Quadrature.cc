#include "LDFE_Quadrature.hh"

namespace // anonymous
{
    int get_number_of_ordinates(dimension,
                                rule)
    {
        vector<int> rules = {4, 16, 64};
    
        switch(rule)
        {
        case 1:
            return rules[0] * 4 * (dimension - 1);
        case 2:
            return rules[1] * 4 * (dimension - 1);
        case 3:
            return rules[2] * 4 * (dimension - 1);
        default:
            AssertMsg(false, "rule not found");
            return -1;
        }
    }
}


LDFE_Quadrature::
LDFE_Quadrature(int dimension,
                int number_of_scattering_moments,
                int rule):
    Angular_Discretization(dimension,
                           number_of_scattering_moments,
                           get_number_of_ordinates(dimension,
                                                   rule)),
    rule_(rule)
{
    initialize_quadrature();
}

void LDFE_Quadrature::
initialize_quadrature()
{
    switch(rule_)
    {
    case 1:
        initialize_1();
        break;
    case 2:
        initialize_2();
        break;
    case 3: 
        initialize_3();
        break;
    }
}

void LDFE_Quadrature::
initialize_1()
{
    vector<double> mu = {};
    vector<double> eta = {};
    vector<double> xi = {};

    ordinates_.resize(number_of_ordinates_ * dimension);
    mu_.resize(number_of_ordinates_);
    eta_.resize(number_of_ordinates_);
    xi_.resize(number_of_ordinates_);
    
    
}

