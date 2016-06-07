#ifndef Line_hh
#define Line_hh

class Line : public Surface
{
public:
    
    Line(vector<double> const &origin,
         vector<double> const &direction);
    
    virtual Relation relation(vector<double> &particle_position) const;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &intersection) const;
    
private:
    
    vector<double> origin_;
    vector<double> direction_;
};

#endif
