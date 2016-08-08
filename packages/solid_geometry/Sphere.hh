#ifndef Sphere_hh
#define Sphere_hh

#include "Surface.hh"

/* 
   Describes a sphere centered at "origin" with radius "radius".
   The formula used is
   
   ||x - x0||^2 = r^2. 
*/
class Sphere : public Surface
{
public:
    
    Sphere(Surface_Type surface_type,
           double radius,
           vector<double> const &origin);
    
    virtual Relation relation(vector<double> const &particle_position) const override;
    virtual double distance(vector<double> const &position) const override;
    virtual Intersection intersection(vector<double> const &particle_position,
                                      vector<double> const &particle_direction,
                                      double &distance,
                                      vector<double> &position) const override;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const override;
private:
    
    double radius_;
    vector<double> origin_;
};

#endif
           
