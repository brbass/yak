#ifndef Line_hh
#define Line_hh

class Line : virtual public Surface
{
public:

    enum class Relation
    {
        POSITIVE,
        NEGATIVE,
        ON
    };

    Line();
    
    virtual Relation relation(double x,
                              double y) const;
    
private:
    
    
    
};

#endif
