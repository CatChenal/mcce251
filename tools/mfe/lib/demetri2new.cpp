#include "mcce.hpp"
#include "math3d.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace Math3D;

int main(int argc, char **argv) {
    Prot prot;
    
    if (argc < 2) {
        cout << "demetri2new in_file_name [out_file_name]" << endl;
        return 0;
    }
    
    prot.load_demetri(argv[1]);
    
    if (argc < 3) {
        prot.write_pdb("step2_out.pdb");
    }
    else {
        prot.write_pdb(argv[2]);
    }
    
    return 0;
}

