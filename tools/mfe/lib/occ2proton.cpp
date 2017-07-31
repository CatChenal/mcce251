#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

int main(int argc, char** argv) {
    Prot prot;
    FILE *fp;
    string msg;
    
    prot.load_conflist("head3.lst");
    prot.load_occtable("fort.38");
    prot.make_sum_crg();
    
    /* writing sum_crg */
    string file_name = "sum_h.out";
    if (!(fp = fopen(file_name.c_str(), "w"))) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    
    if (prot.param.titra_type == PH) {
        fprintf(fp, " ph       ");
        for (int i=0; i<prot.titra_point.size(); i++)
            fprintf(fp, " %5.1f", prot.titra_point[i]);
    }
    else {
        fprintf(fp, " eh       ");
        for (int i=0; i<prot.titra_point.size(); i++)
            fprintf(fp, " %5.0f", prot.titra_point[i]);
    }
    fprintf(fp, "\n");
    
    for (int i_res=0;i_res<prot.res.size();++i_res) {
        int i_conf;
        for (i_conf=1;i_conf<prot.res[i_res].conf.size();++i_conf) {
            if (fabs(prot.res[i_res].conf[i_conf].H) > 1e-4) break;
        }
        if (i_conf >= prot.res[i_res].conf.size()) continue;
        
        fprintf(fp, "%s ", prot.res[i_res].resID().c_str());
        fprintf(fp, "%s\n", vec2string(prot.res[i_res].sum_H, "%5.2f").c_str());
    }
    
    fprintf(fp, "Tot Protn ");
    fprintf(fp, "%s\n", vec2string(prot.sum_H, "%5.1f").c_str());
    fclose(fp);
    return 0;
}

