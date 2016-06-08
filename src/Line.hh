#ifndef Line_hh
#define Line_hh

class Line : public Surface
{
public:
    
    Line(Surface_Type surface_type,
         vector<double> const &origin,
         vector<double> const &direction,
         vector<double> const &normal);
    
    virtual Relation relation(vector<double> &particle_position) const;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal) const;
    
private:
    
    vector<double> origin_;
    vector<double> direction_;
    vector<double> normal_;
};

#endif
