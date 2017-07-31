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
        cout << "   vdw whole_pdb small_pdb\n";
        exit(-1);
    }
    
    prot.load_pdb(argv[1]);
    prot_small.load_pdb(argv[2]);
    
    ofstream conf_vdw;
    conf_vdw.open("conf_vdw.out", ios::out);

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
                                if (atom1->atomName[1] == 'C' && atom2->atomName[1] == 'C')
                                    e += 2516582.400/pow(d,12) - 1228.800000/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'C' && atom2->atomName[1] == 'N')
                                    e += 1198066.249/pow(d,12) - 861.634784/pow(d,6);
                                else if (atom1->atomName[1] == 'N' && atom2->atomName[1] == 'C')
                                    e += 1198066.249/pow(d,12) - 861.634784/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'C' && atom2->atomName[1] == 'O')
                                    e += 820711.722/pow(d,12) - 754.059521/pow(d,6);
                                else if (atom1->atomName[1] == 'O' && atom2->atomName[1] == 'C')
                                    e += 820711.722/pow(d,12) - 754.059521/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'C' && atom2->atomName[1] == 'S')
                                    e += 2905899.052/pow(d,12) - 1418.896022/pow(d,6);
                                else if (atom1->atomName[1] == 'S' && atom2->atomName[1] == 'C')
                                    e += 2905899.052/pow(d,12) - 1418.896022/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'C' && atom2->atomName[1] == 'H')
                                    e += 29108.222/pow(d,12) - 79.857949/pow(d,6);
                                else if (atom1->atomName[1] == 'H' && atom2->atomName[1] == 'C')
                                    e += 29108.222/pow(d,12) - 79.857949/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'N' && atom2->atomName[1] == 'N')
                                    e += 540675.281/pow(d,12) - 588.245000/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'N' && atom2->atomName[1] == 'O')
                                    e += 357365.541/pow(d,12) - 505.677729/pow(d,6);
                                else if (atom1->atomName[1] == 'O' && atom2->atomName[1] == 'N')
                                    e += 357365.541/pow(d,12) - 505.677729/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'N' && atom2->atomName[1] == 'S')
                                    e += 1383407.742/pow(d,12) - 994.930149/pow(d,6);
                                else if (atom1->atomName[1] == 'S' && atom2->atomName[1] == 'N')
                                    e += 1383407.742/pow(d,12) - 994.930149/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'N' && atom2->atomName[1] == 'H')
                                    e += 10581.989/pow(d,12) - 48.932922/pow(d,6);
                                else if (atom1->atomName[1] == 'H' && atom2->atomName[1] == 'N')
                                    e += 10581.989/pow(d,12) - 48.932922/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'O' && atom2->atomName[1] == 'O')
                                    e += 230584.301/pow(d,12) - 429.496730/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'O' && atom2->atomName[1] == 'S')
                                    e += 947676.268/pow(d,12) - 870.712934/pow(d,6);
                                else if (atom1->atomName[1] == 'S' && atom2->atomName[1] == 'O')
                                    e += 947676.268/pow(d,12) - 870.712934/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'O' && atom2->atomName[1] == 'H')
                                    e += 6035.457/pow(d,12) - 39.075098/pow(d,6);
                                else if (atom1->atomName[1] == 'H' && atom2->atomName[1] == 'O')
                                    e += 6035.457/pow(d,12) - 39.075098/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'S' && atom2->atomName[1] == 'S')
                                    e += 3355443.200/pow(d,12) - 1638.400000/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'S' && atom2->atomName[1] == 'H')
                                    e += 33611.280/pow(d,12) - 92.212017/pow(d,6);
                                else if (atom1->atomName[1] == 'H' && atom2->atomName[1] == 'S')
                                    e += 33611.280/pow(d,12) - 92.212017/pow(d,6);
                                
                                else if (atom1->atomName[1] == 'H' && atom2->atomName[1] == 'H')
                                    e += 81.920/pow(d,12) - 2.560000/pow(d,6);
                                
                                else
                                    e += 2905899.052/pow(d,12) - 1418.896022/pow(d,6);
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
                conf_vdw << prot.res[i_res].conf[i_conf].uniqID << " OVERF\n";
            }
            else {
                char e_str[80];
                sprintf(e_str,"%6.2f",conf_e);
                conf_vdw << prot.res[i_res].conf[i_conf].uniqID << " " << e_str << endl;
            }
        }
    }
}

/* 
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H1$
VDWAMBER C12   C-C  2516582.400     # reqm = 4.00, eij = -0.150 kcal/mol
VDWAMBER C6    C-C  1228.800000
VDWAMBER C12   C-N  1198066.249     # reqm = 3.75, eij = -0.155 kcal/mol
VDWAMBER C6    C-N   861.634784
VDWAMBER C12   C-O   820711.722     # reqm = 3.60, eij = -0.173 kcal/mol
VDWAMBER C6    C-O   754.059521
VDWAMBER C12   C-S  2905899.052     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    C-S  1418.896022
VDWAMBER C12   C-H    29108.222     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    C-H    79.857949

VDWAMBER C12   N-C  1198066.249     # reqm = 3.75, eij = -0.155 kcal/mol
VDWAMBER C6    N-C   861.634784
VDWAMBER C12   N-N   540675.281     # reqm = 3.50, eij = -0.160 kcal/mol
VDWAMBER C6    N-N   588.245000
VDWAMBER C12   N-O   357365.541     # reqm = 3.35, eij = -0.179 kcal/mol
VDWAMBER C6    N-O   505.677729
VDWAMBER C12   N-S  1383407.742     # reqm = 3.75, eij = -0.179 kcal/mol
VDWAMBER C6    N-S   994.930149
VDWAMBER C12   N-H    10581.989     # reqm = 2.75, eij = -0.057 kcal/mol
VDWAMBER C6    N-H    48.932922

VDWAMBER C12   O-C   820711.722     # reqm = 3.60, eij = -0.173 kcal/mol
VDWAMBER C6    O-C   754.059521
VDWAMBER C12   O-N   357365.541     # reqm = 3.35, eij = -0.179 kcal/mol
VDWAMBER C6    O-N   505.677729
VDWAMBER C12   O-O   230584.301     # reqm = 3.20, eij = -0.200 kcal/mol
VDWAMBER C6    O-O   429.496730
VDWAMBER C12   O-S   947676.268     # reqm = 3.60, eij = -0.200 kcal/mol
VDWAMBER C6    O-S   870.712934
VDWAMBER C12   O-H     6035.457     # reqm = 2.60, eij = -0.063 kcal/mol
VDWAMBER C6    O-H    39.075098

VDWAMBER C12   S-C  2905899.052     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    S-C  1418.896022
VDWAMBER C12   S-N  1383407.742     # reqm = 3.75, eij = -0.179 kcal/mol
VDWAMBER C6    S-N   994.930149
VDWAMBER C12   S-O   947676.268     # reqm = 3.60, eij = -0.200 kcal/mol
VDWAMBER C6    S-O   870.712934
VDWAMBER C12   S-S  3355443.200     # reqm = 4.00, eij = -0.200 kcal/mol
VDWAMBER C6    S-S  1638.400000
VDWAMBER C12   S-H    33611.280     # reqm = 4.00, eij = -0.200 kcal/mol
VDWAMBER C6    S-H    92.212017

VDWAMBER C12   H-C    29108.222     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    H-C    79.857949
VDWAMBER C12   H-N    10581.989     # reqm = 2.75, eij = -0.057 kcal/mol
VDWAMBER C6    H-N    48.932922
VDWAMBER C12   H-O     6035.457     # reqm = 2.60, eij = -0.063 kcal/mol
VDWAMBER C6    H-O    39.075098
VDWAMBER C12   H-S    33611.280     # reqm = 4.00, eij = -0.200 kcal/mol
VDWAMBER C6    H-S    92.212017
VDWAMBER C12   H-H       81.920     # reqm = 2.00, eij = -0.020 kcal/mol
VDWAMBER C6    H-H     2.560000

# default is the same as S-C
VDWAMBER C12   X-X  2905899.052     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    X-X  1418.896022
 */
