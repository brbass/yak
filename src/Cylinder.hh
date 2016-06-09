#ifndef Cylinder_hh
#define Cylinder_hh

#include "Surface.hh"

class Cylinder : public Surface
{
public:
    
    Cylinder(Surface_Type surface_type,
             double radius,
             vector<double> const &origin,
             vector<double> const &direction);
    
    virtual Relation relation(vector<double> const &particle_position) const;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const;
    
protected:
    
    double radius_;
    vector<double> origin_;
    vector<double> direction_;
};

#endif
