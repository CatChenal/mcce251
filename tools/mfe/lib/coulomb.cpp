#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "mcce.hpp"
#include "math3d.hpp"

using namespace std;

int main(int argc, char **argv) {
    Prot prot, prot_small;
    
    if (argc < 3) {
        cout << "Usage: \n";
        cout << "   coulumb whole_pdb small_pdb\n";
        exit(-1);
    }
    
    prot.load_pdb(argv[1]);
    prot_small.load_pdb(argv[2]);
    
    ofstream conf_coulomb;
    conf_coulomb.open("conf_coulomb.out", ios::out);

    for (int i_res=0; i_res<prot.res.size(); ++i_res) {
        for (int i_conf=0; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
            float conf_e = 0;
            for (int i_atom=0; i_atom<prot.res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                Atom* atom1 = &prot.res[i_res].conf[i_conf].atom[i_atom];
                
                float e = 0;
                
                for (int j_res=0; j_res<prot_small.res.size(); ++j_res) {
                    for (int j_conf=0; j_conf<prot_small.res[j_res].conf.size(); ++j_conf) {
                        for (int j_atom=0; j_atom<prot_small.res[j_res].conf[j_conf].atom.size(); ++j_atom) {
                            Atom* atom2 = &prot_small.res[j_res].conf[j_conf].atom[j_atom];
                            
                            float d = (atom1->r - atom2->r).Length3();
                            
                            if (d>0.5) {
                                e += 331.5*atom1->crg*atom2->crg/(4. * d);
                            }
                            else {
                                e = 1000;
                                break;
                            }
                        }
                        if (e>999.99) break;
                    }
                    if (e>999.99) break;
                }
                
                string pdbstr;
                atom1->write_pdbstr(pdbstr);
                if (e>999.99) {
                    pdbstr.replace(6,5,"OVERF");
                    conf_e = 1000.;
                    break;
                }
                else {
                    char e_str[80];
                    sprintf(e_str,"%6.2f",e);
                    pdbstr.replace(5,6,e_str);
                    conf_e += e;
                }
                cout << pdbstr << endl;
            }
            
            if (conf_e>999.99) {
                conf_coulomb << prot.res[i_res].conf[i_conf].uniqID << " OVERF\n";
            }
            else {
                char e_str[80];
                sprintf(e_str,"%6.2f",conf_e);
                conf_coulomb << prot.res[i_res].conf[i_conf].uniqID << " " << e_str << endl;
            }
        }
    }
}
