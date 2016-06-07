#ifndef Axis_Line_hh
#define Axis_Line_hh

class Axis_Line : public Line
{
public:

    Axis_Line(int dimension,
              double value);
    
    virtual Relation relation(double x,
                              double y) const;
    
private:

    int dimension_;
    double value_;
};

#endif
