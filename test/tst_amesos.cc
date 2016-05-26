#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Amesos.h>

using namespace std;

int main()
{
    Amesos factory;

    vector<string> solver_types = {"Btf",
                                   "Dscpack",
                                   "Klu",
                                   "Lapack",
                                   "Merikos",
                                   "Mumps",
                                   "Paraklete",
                                   "Pardiso",
                                   "Scalapack",
                                   "Superlu",
                                   "Superludist",
                                   "Taucs",
                                   "Umfpack"};
    
    for (string st : solver_types)
    {
        cout << setw(16) << st;
        cout << setw(16) << factory.Query(st.c_str());
        cout << endl;
    }
}
