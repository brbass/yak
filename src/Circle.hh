#ifndef Circle_hh
#define Circle_hh

class Circle : public Surface
{
public:

    Circle(double radius,
           vector<double> const &origin);

    virtual Relation relation(vector<double> &particle_position) const;
    virtual bool intersection(vector<double> const &particle_position,
                              vector<double> const &particle_direction,
                              double &distance,
                              vector<double> &intersection) const;
    
private:

    double radius_;
    vector<double> origin_;
};

#endif
           
