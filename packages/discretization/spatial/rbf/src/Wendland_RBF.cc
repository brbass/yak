#include "Wendland_RBF.hh"

#include <cmath>
#include <string>

#include "Check.hh"

using std::to_string;

Wendland_RBF::
Wendland_RBF(int number_of_dimensions,
             int order,
             vector<double> const &position,
             vector<double> const &shape_parameter):
    RBF(number_of_dimensions,
        position,
        shape_parameter),
    order_(order)
{
}

double Wendland_RBF::
basis(vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double dist = sqrt(dist2);

    switch(order_)
    {
    case 0:
        if(dist<1)
        {
            return pow(1-dist, 2);
        }
        return 0;
    case 1:
        if (dist<1)
        {
            return pow(1-dist, 4)*(4*dist+1);
        }
        return 0;
    case 2:
        if (dist<1)
        {
            return pow(1-dist, 6)*(35*dist*dist+18*dist+3);
        }
        return 0;
    case 3:
        if (dist<1)
        {
            return pow(1-dist, 8)*(32*dist*dist*dist+18*dist*dist+8*dist+1);
        }
        return 0;
    default:
        AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        
        return 0;
    }
}

double Wendland_RBF::
dbasis(int dim,
       vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double dist = sqrt(dist2);
    double kr2 = shape_[dim] * shape_[dim] * (r[dim] - position_[dim]);
    
    switch(order_)
    {
    case 0:
        if (dist==0.)
        {
            return -2 * shape(dim);
        }
        else if(dist<1)
        {
            return kr2*2*(1-1/dist);
        }
        return 0;
    case 1:
        if (dist<1)
        {
            return kr2*20*pow(dist-1, 3);
        }
        return 0;
    case 2:
        if (dist<1)
        {
            return kr2*56*pow(dist-1, 5)*(5*dist+1);
        }
        return 0;
    case 3:
        if (dist<1)
        {
            return kr2*22*pow(dist-1, 7)*(1+dist*(16*dist+7));
        }
        return 0;
    default:
        AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        
        return 0;
    }
}

double Wendland_RBF::
ddbasis(int dim,
        vector<double> const &r) const
{
    double dist2 = get_distance_squared(r);
    double dist = sqrt(dist2);
    double distdim = r[dim] - position_[dim];
    double k2r2 = shape_[dim] * shape_[dim] * distdim * distdim;

    switch(order_)
    {
    case 0:
        if (dist==0.)
        {
            return -2 * shape_[dim] * shape_[dim];
        }
        else if(dist<1)
        {
            return 2 * shape_[dim] * shape_[dim] * (dist2 * (dist - 1) + k2r2) / (dist2 * dist);
        }
        return 0;
    case 1:
        if (dist==0.)
        {
            return -20 * shape_[dim] * shape_[dim];
        }
        else if (dist<1)
        {
            return 20 * shape_[dim] * shape_[dim] * pow(dist - 1, 2)* (dist*(dist-1) + 3 * k2r2) / dist;
        }
        return 0;
    case 2:
        if (dist<1)
        {
            return 56./3.*shape_[dim] * shape_[dim] * pow(dist-1, 4)*(-1 + dist*(-4+5*dist)+30*k2r2);
        }
        return 0;
    case 3:
        if (dist<1)
        {
            return 22 * shape_[dim] * shape_[dim]* pow(dist - 1, 6)*((dist-1)*(1+dist*(7+16*dist)) + 24*k2r2*(1+6*dist));
        }
        return 0;
    default:
        AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        
        return 0;
    }
}
