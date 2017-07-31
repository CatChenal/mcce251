#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "mcce.hpp"
#include "math3d.hpp"

using namespace std;

int main(int argc, char **argv) {
    Prot prot;
    
    if (argc < 2) {
        cout << "Usage: \n";
        cout << "   fix_boundary file.pdb\n";
        exit(-1);
    }

    prot.load_pdb(argv[1]);
    prot.index_on();
    int shift_flag = 1;
    while (shift_flag) {
        shift_flag = 0;
        for (int i_res=0; i_res<prot.res.size(); ++i_res) {
            printf("%d\n",i_res);
            for (int i_conf=0; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
                for (int i_atom=0; i_atom<prot.res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                    Atom* atom1 = &prot.res[i_res].conf[i_conf].atom[i_atom];
                    
                    for (int j_res=i_res; j_res<prot.res.size(); ++j_res) {
                        for (int j_conf=0; j_conf<prot.res[j_res].conf.size(); ++j_conf) {
                            for (int j_atom=0; j_atom<prot.res[j_res].conf[j_conf].atom.size(); ++j_atom) {
                                Atom* atom2 = &prot.res[j_res].conf[j_conf].atom[j_atom];
                                if (i_res==j_res && i_conf==j_conf && i_atom==j_atom)continue;
                                
                                if ( (atom1->r - atom2->r).Length3() > 0.002) continue;
                                if (fabs(atom1->rad -atom2->rad) < 0.01) continue;
                                
                                shift_flag =1;
                                if (atom1->rad < atom2->rad) atom1->r(2) += 0.002;
                                else atom2->r(2) += 0.002;
                                
                                string pdbstr;
                                atom1->write_pdbstr(pdbstr);
                                printf("%s\n", pdbstr.c_str());
                                atom2->write_pdbstr(pdbstr);
                                printf("%s\n\n", pdbstr.c_str());
                            }
                        }
                    }
                }
            }
        }
    }
    
    prot.write_pdb("new.pdb");
}
