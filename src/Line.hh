#ifndef Line_hh
#define Line_hh

#include "Surface.hh"

/* 
   Describes an infinite line that goes through the point "origin"
   with normal vector "normal". The equation is
   
   \Omega \cdot (x - x0) = 0.
*/
class Line : public Surface
{
public:
    
    Line(Surface_Type surface_type,
         vector<double> const &origin,
         vector<double> const &normal);
    
    virtual Relation relation(vector<double> const &particle_position) const;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const;
    
protected:
    
    vector<double> origin_;
    vector<double> normal_;
};

#endif
