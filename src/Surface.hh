#ifndef Surface_hh
#define Surface_hh

class Surface
{
public:

    enum class Relation
    {
        POSITIVE,
        NEGATIVE,
        ON
    }
    
    Surface();

private:

    virtual Relation relation(double x) const;
    virtual Relation relation(double x,
                              double y) const;
    virtual Relation relation(double x,
                              double y,
                              double z) const;

}

#endif
