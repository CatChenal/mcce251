#include "mcce.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int db_open();

int main(int argc, char *argv[]) {
    ifstream inpdb;
    ofstream outpdb;
    Prot prot;
    
    try {
        prot.load_pdb(argv[1]);
    }
    catch(string msg) {
        cout << "exiting!" << endl;
        exit(-1);
    }
    
    cout << prot.count_conf() << " confs " << prot.count_atom() << " atoms" << endl;
    
    
    try {
        prot.write_pdb(argv[2]);
    }
    catch(string msg) {
        cout << "exiting!" << endl;
        exit(-1);
    }
    
    int i_res, i_conf, i_atom;
    for (i_res = 0; i_res < prot.res.size(); ++i_res) {
        for (i_conf = 0; i_conf < prot.res[i_res].conf.size(); ++i_conf) {
            for (i_atom = 0; i_atom < prot.res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                if (prot.res[i_res].conf[i_conf].atom[i_atom].off()) continue;
                prot.res[i_res].conf[i_conf].atom[i_atom].r.Normalize3();
            }
        }
    }
    return 0;
}

