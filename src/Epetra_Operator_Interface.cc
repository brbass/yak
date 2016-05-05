#include "Check.hh"
#include "Epetra_Operator_Interface.hh"

Epetra_Operator_Interface::
Epetra_Operator_Interface(shared_ptr<Epetra_Comm> const &comm,
                          shared_ptr<Epetra_Map> const &map,
                          shared_ptr<Vector_Operator> const &oper):
    comm_(comm),
    map_(map),
    oper_(oper)
{
    Check(oper_->square());
}

Epetra_Operator_Interface::
~Epetra_Operator_Interface()
{
}

int Epetra_Operator_Interface::
Apply(Epetra_MultiVector const &X,
      Epetra_MultiVector &Y) const
{
    Assert(X.NumVectors() == Y.NumVectors());
    Assert(X.NumVectors() == 1);
    int const size = oper_->row_size();
    
    vector<double> x(size);
    X.ExtractCopy(&x[0], X.Stride());
    
    (*oper_)(x);

    for (int i = 0; i < size; ++i)
    {
        Y.ReplaceMyValue(i, 0, x[i]);
    }
}
