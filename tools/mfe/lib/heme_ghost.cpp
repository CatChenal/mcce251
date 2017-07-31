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
    
    if (argc < 3) {
        cout << "heme_ghost in_file_name distance" << endl;
        return 0;
    }
    
    prot.load_pdb(argv[1]);
    for (int i_res=0; i_res< prot.res.size(); ++i_res) {
        if (prot.res[i_res].resName == "HA3"
        || prot.res[i_res].resName == "HMB" ) {
            for (int i_conf=1; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
                //if (prot.res[i_res].conf[i_conf].confName.substr(3,2) == "+W") {
                    Vector4<float> r_a, r_b, r_c, r_d, r_ac, r_bd, r_fe, r_his, r_polar;
                    int assigned = 0;
                    for (int i_atom=0; i_atom<prot.res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                        string atomname = strip(prot.res[i_res].conf[i_conf].atom[i_atom].atomName);
                        if ( atomname == "N A") {
                            prot.res[i_res].conf[i_conf].atom[i_atom].crg = 1.0;
                            r_a = prot.res[i_res].conf[i_conf].atom[i_atom].r;
                            assigned++;
                        }
                        if ( atomname == "N B") {
                            prot.res[i_res].conf[i_conf].atom[i_atom].crg = 1.0;
                            r_b = prot.res[i_res].conf[i_conf].atom[i_atom].r;
                            assigned++;
                        }
                        if ( atomname == "N C") {
                            prot.res[i_res].conf[i_conf].atom[i_atom].crg = 1.0;
                            r_c = prot.res[i_res].conf[i_conf].atom[i_atom].r;
                            assigned++;
                        }
                        if ( atomname == "N D") {
                            prot.res[i_res].conf[i_conf].atom[i_atom].crg = 1.0;
                            r_d = prot.res[i_res].conf[i_conf].atom[i_atom].r;
                            assigned++;
                        }
                        if ( atomname == "FE") {
                            prot.res[i_res].conf[i_conf].atom[i_atom].crg = -1.0;
                            r_fe = prot.res[i_res].conf[i_conf].atom[i_atom].r;
                            assigned++;
                        }
                        if ( atomname == "NE2") {
                            r_his = prot.res[i_res].conf[i_conf].atom[i_atom].r;
                            assigned++;
                        }
                    }
                    if (assigned<6) {
                        cout << " Only " << assigned << " out of 6 atoms are found" << endl;
                        continue;
                    }
                    
                    r_ac = (r_c - r_a);
                    r_bd = (r_d - r_b);
                    r_polar = (CrossProduct(r_ac, r_bd));
                    r_polar.Normalize3();
                    if (r_polar.DotProduct3(r_fe - r_his) < 0.) {
                        Vector4<float> zero(0,0,0);
                        r_polar = zero - r_polar;
                    }
                    
                    r_polar *= atof(argv[2]);
                    
                    Atom new_atom;
                    new_atom = prot.res[i_res].conf[i_conf].atom[0];
                    
                    new_atom.rad = 1.00;
                    new_atom.crg = -0.25;
                    new_atom.atomName = "BQA1";
                    new_atom.r = r_a + r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                    
                    new_atom.atomName = "BQB1";
                    new_atom.r = r_b + r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                    
                    new_atom.atomName = "BQC1";
                    new_atom.r = r_c + r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                    
                    new_atom.atomName = "BQD1";
                    new_atom.r = r_d + r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
    
                    new_atom.crg = -0.25;
                    new_atom.atomName = "BQA2";
                    new_atom.r = r_a - r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                    
                    new_atom.atomName = "BQB2";
                    new_atom.r = r_b - r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                    
                    new_atom.atomName = "BQC2";
                    new_atom.r = r_c - r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                    
                    new_atom.atomName = "BQD2";
                    new_atom.r = r_d - r_polar;
                    prot.res[i_res].conf[i_conf].ins_atom(new_atom);
                //}
            }
        }
    }
    prot.write_pdb("Heme_ghost.pdb");
    return 0;
}

