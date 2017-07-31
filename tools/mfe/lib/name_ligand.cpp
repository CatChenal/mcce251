#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

int main() {
    Prot prot;
    
    prot.load_pdb("prot.pdb");
    
    for (int i_res=0; i_res <prot.res.size(); ++i_res) {
        if (prot.res[i_res].resName != "HEM") continue;
        int i_atom;
        for (i_atom=0;i_atom<prot.res[i_res].conf[1].atom.size(); ++i_atom) {
            if (prot.res[i_res].conf[1].atom[i_atom].atomName == "FE  ") break;
        }
        if (i_atom>=prot.res[i_res].conf[1].atom.size()){
            cerr << "Can't find Fe\n";
            return -1;
        }

        for (int j_res=0; j_res <prot.res.size(); ++j_res) {
            if (prot.res[j_res].resName != "HIS") continue;
            int j_atom;
            for (j_atom=0; j_atom<prot.res[j_res].conf[1].atom.size(); ++j_atom) {
                if (prot.res[j_res].conf[1].atom[j_atom].atomName == " NE2") break;
            }
            if (i_atom>=prot.res[i_res].conf[1].atom.size()){
                continue;
            }
            
            if ((prot.res[i_res].conf[1].atom[i_atom].r - prot.res[j_res].conf[1].atom[j_atom].r).Length3() < 3.0) {
                FILE *fp;
                fp = fopen("name.txt","a");
                if (prot.res[i_res].chainID == ' ') prot.res[i_res].chainID = '_';
                if (prot.res[j_res].chainID == ' ') prot.res[j_res].chainID = '_';
                
                fprintf(fp," CA *HIS*%c%4d  *****BKB******\n",prot.res[j_res].chainID, prot.res[j_res].resSeq);
                fprintf(fp," N  *HIS*%c%4d  *****BKB******\n",prot.res[j_res].chainID, prot.res[j_res].resSeq);
                fprintf(fp," C  *HIS*%c%4d  *****BKB******\n",prot.res[j_res].chainID, prot.res[j_res].resSeq);
                fprintf(fp," O  *HIS*%c%4d  *****BKB******\n",prot.res[j_res].chainID, prot.res[j_res].resSeq);
                
                fprintf(fp, "*****HIS*%c%4d  *****HAN*%c%04d\n",
                prot.res[j_res].chainID, prot.res[j_res].resSeq,
                prot.res[i_res].chainID, prot.res[i_res].resSeq);

                fclose(fp);
            }
        }
    }
    return 0;
}
