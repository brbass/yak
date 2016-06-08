#ifndef Line_hh
#define Line_hh

class Line : public Surface
{
public:
    
    Line(bool boundary,
         double reflection_probability,
         vector<double> const &origin,
         vector<double> const &direction);
    
    virtual Relation relation(vector<double> &particle_position) const;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &position) const;
    
private:
    
    vector<double> origin_;
    vector<double> direction_;
};

#endif
