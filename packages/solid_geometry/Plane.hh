#ifndef Plane_hh
#define Plane_hh

#include "Surface.hh"

/* 
   Describes an infinite plane that goes through the point "origin" 
   and has surface normal "normal". The formula used is
   
   \Omega \cdot (x - x0) = 0
*/
class Plane : public Surface
{
public:
    
    Plane(Surface_Type surface_type,
          vector<double> const &origin,
          vector<double> const &normal);
    
    virtual Relation relation(vector<double> const &particle_position) const override;
    virtual Intersection intersection(vector<double> const &particle_position,
                                      vector<double> const &particle_direction,
                                      double &distance,
                                      vector<double> &position) const override;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const override;
    
protected:
    
    vector<double> origin_;
    vector<double> normal_;
};

#endif
