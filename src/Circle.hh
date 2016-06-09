#ifndef Circle_hh
#define Circle_hh

#include "Surface.hh"

/*
  Describes a circle with radius "radius" centered at the point "origin."
  The formula used is
  
  ||x - x0||^2 = r^2.
*/
class Circle : public Surface
{
public:
    
    Circle(Surface_Type surface_type,
           double radius,
           vector<double> const &origin);
    
    virtual Relation relation(vector<double> const &particle_position) const;
    virtual Intersection intersection(vector<double> const &particle_position,
                                      vector<double> const &particle_direction,
                                      double &distance,
                                      vector<double> &position) const;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const;
private:
    
    double radius_;
    vector<double> origin_;
};

#endif
           
