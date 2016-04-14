#ifndef Epetra_Operator_Interface_hh
#define Epetra_Operator_Interface_hh

#include <memory>

#ifdef EPETRA_MPI
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include "Vector_Operator.hh"

using std::shared_ptr;

class Epetra_Operator_Interface: public Epetra_Operator
{
public:

    Epetra_Operator_Interface(shared_ptr<Epetra_Comm> const &comm,
                              shared_ptr<Epetra_Map> const &map,
                              shared_ptr<Vector_Operator> const &oper);

    ~Epetra_Operator_Interface();

    virtual int SetUseTranspose(bool UseTranspose)
    {
        return -1;
    }

    virtual int Apply(Epetra_MultiVector const &X,
                      Epetra_MultiVector &Y) const;

    virtual int ApplyInverse(Epetra_MultiVector const &X,
                             Epetra_MultiVector &Y) const
    {
        return 1;
    }

    virtual double NormInf() const
    {
        return 0.;
    }

    virtual const char *Label() const
    {
        return "Epetra_Interface";
    }

    virtual bool UseTranspose() const
    {
        return false;
    }

    virtual bool HasNormInf() const
    {
        return false;
    }

    virtual const Epetra_Comm &Comm() const
    {
        return *comm_;
    }
    
    virtual const Epetra_Map &OperatorDomainMap() const
    {
        return *map_;
    }

    virtual const Epetra_Map &OperatorRangeMap() const
    {
        return *map_;
    }
    

private:
    
    shared_ptr<Epetra_Comm> comm_;
    shared_ptr<Epetra_Map> map_;
    shared_ptr<Vector_Operator> oper_;
};

#endif
