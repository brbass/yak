#include <iostream>
#include <string>

#include "Driver.hh"

using namespace std;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cerr << "usage: yak [input.xml]" << endl;
        return 1;
    }
    
    string filename = argv[1];
    
    Driver driver(filename);
}
