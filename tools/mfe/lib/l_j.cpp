#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "mcce.hpp"
#include "math3d.hpp"

using namespace std;

double vdw(Conf &conf1, Conf &conf2);
double vdw_amber(Conf &conf1, Conf &conf2);
double coulomb(Conf &conf1, Conf &conf2);

int main(int argc, char **argv) {
    Prot prot;
    Atom *atom1, *atom2;
    
    if (argc < 2) {
        cout << "Usage: \n";
        cout << "   " << argv[0] << " pdb\n";
        exit(-1);
    }
    
    prot.load_pdb(argv[1]);
    
    for (int k_atom = 0; k_atom < prot.res[0].conf[1].atom.size(); ++k_atom) {
        if (prot.res[0].conf[1].atom[k_atom].atomName == " OD2") {
            atom1 = &prot.res[0].conf[1].atom[k_atom];
            break;
        }
    }
    for (int k_atom = 0; k_atom < prot.res[1].conf[1].atom.size(); ++k_atom) {
        if (prot.res[1].conf[1].atom[k_atom].atomName == " OH ") {
            atom2 = &prot.res[1].conf[1].atom[k_atom];
            break;
        }
    }
    
    Vect translation;
    translation = (atom2->r-atom1->r).Normalize3();
    for (int k_atom = 0; k_atom < prot.res[1].conf[1].atom.size(); ++k_atom) {
        prot.res[1].conf[1].atom[k_atom].r -= translation;
    }
    
    translation *= 0.1;
    for (int i_step = 0; i_step <31; ++i_step) {
        for (int k_atom = 0; k_atom < prot.res[1].conf[1].atom.size(); ++k_atom) {
            prot.res[1].conf[1].atom[k_atom].r += translation;
        }
        printf("%8.3f %11.3f %11.3f %11.3f\n",
        (atom1->r - atom2->r).Length3(),
        vdw_amber(prot.res[0].conf[1],prot.res[1].conf[1]),
        coulomb(prot.res[0].conf[1],prot.res[1].conf[1]),
        vdw_amber(prot.res[0].conf[1],prot.res[1].conf[1])+coulomb(prot.res[0].conf[1],prot.res[1].conf[1]));
    }
    return 0;
}

double vdw(Conf &conf1, Conf &conf2)
{
    double e = 0;
    for (int i_atom=0; i_atom < conf1.atom.size(); ++i_atom) {
        Atom* atom1 = &conf1.atom[i_atom];
        
        for (int j_atom=0; j_atom<conf2.atom.size(); ++j_atom) {
            Atom* atom2 = &conf2.atom[j_atom];
            
            double d = (atom1->r - atom2->r).Length3();
            
            if (d<0.5) d=0.5;
            
            /* amber autodock parameter */
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
    }
    
    return e;
}

double coulomb(Conf &conf1, Conf &conf2)
{
    double e = 0;
    for (int i_atom=0; i_atom<conf1.atom.size(); ++i_atom) {
        Atom* atom1 = &conf1.atom[i_atom];
        
        for (int j_atom=0; j_atom<conf2.atom.size(); ++j_atom) {
            Atom* atom2 = &conf2.atom[j_atom];
            
            double d = (atom1->r - atom2->r).Length3();
            
            e += 331.5*atom1->crg*atom2->crg/d;
        }
    }
    
    return e;
}

double vdw_amber(Conf &conf1, Conf &conf2)
{
    double e = 0;
    for (int i_atom=0; i_atom < conf1.atom.size(); ++i_atom) {
        Atom* atom1 = &conf1.atom[i_atom];
        
        for (int j_atom=0; j_atom<conf2.atom.size(); ++j_atom) {
            Atom* atom2 = &conf2.atom[j_atom];
            
            double d = (atom1->r - atom2->r).Length3();
            
            if (d<0.5) d=0.5;
            
            /* amber autodock parameter */
            if (atom1->atomName == " OD1" && atom2->atomName == " HH ")
                e += 0.;
            else if (atom1->atomName == " OD2" && atom2->atomName == " HH ")
                e += 0.;

            else if (atom1->atomName == " OD1" && atom2->atomName == " OH ")
                e += 471003.2866/pow(d,12) - 629.3007104/pow(d,6);
            else if (atom1->atomName == " OD2" && atom2->atomName == " OH ")
                e += 471003.2866/pow(d,12) - 629.3007104/pow(d,6);
            
            else if (atom1->atomName[1] == 'C' && atom2->atomName[1] == 'C')
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
    }
    
    return e;
}

