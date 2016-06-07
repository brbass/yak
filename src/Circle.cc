#include "Circle.hh"

#include <cmath>

namespace // anonymous
{
    double dot(vector<double> const &x,
               vector<double> const &y)
    {
        return x[0] * y[0] + x[1] * y[1];
    }
    
    double cross(vector<double> const &x,
                 vector<double> const &y)
    {
        return x[0] * y[1] - x[1] * y[0];
    }
}

Circle::
Circle(double radius,
       vector<double> const &origin):
    radius_(radius),
    origin_(origin)
{
}

Relation Circle::
relation(vector<double> &particle_position) const
{
    double x = particle_position[0] - origin_[0];
    double y = particle_position[1] - origin_[1];
    
    double r = sqrt(x * x + y * y);
    
    if (r < radius_)
    {
        return Relation::INSIDE;
    }
    else if (r > radius_)
    {
        return Relation::OUTSIDE;
    }
    else // r == radius_
    {
        return Relation::EQUAL;
    }
}

bool Circle::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    vector<double> center_distance(dimension_);
    
    for (int i = 0; i < dimension_; ++i)
    {
        center_distance[i] = origin_[i] - particle_position[i];
    }

    double cross_val = cross(particle_direction, center_distance);
    double k1 = radius_ * radius_ - cross_val * cross_val;

    if (k1 <= 0)
    {
        return false;
    }

    double k2 = sqrt(k1);
    double k3 = dot(center_distance, particle_direction);
    
    double t1 = k3 + k2;
    double t2 = k3 - k2;
    double t;
    
    if (t2 > 0)
    {
        t = t2;
    }
    else if (t1 > 0)
    {
        t = t1;
    }
    else
    {
        return false;
    }

    distance = t;
    
    position.resize(dimension_);
    
    for (int i = 0; i < dimension_; ++i)
    {
        position[i] = particle_position[i] + particle_direction[i] * t;
    }
    
    return true;
}
