#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

int main() {
    Prot prot;
    
    prot.load_conflist("head3.lst");
    prot.load_occtable("fort.38");
    prot.index_on();
    
    vector<float> sum_occ(prot.titra_point.size());
    
    for (int i_res=0; i_res<prot.res.size(); ++i_res) {
        if (prot.res[i_res].resName != "HA3") continue;
        for (int i_conf=1; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
            if (prot.res[i_res].conf[i_conf].uniqID.substr(0,4) == "HA3+") {
                sum_occ = sum_2vec(sum_occ, prot.res[i_res].conf[i_conf].occ_table);
            }
        }
        cout << " eh        " << vec2string(prot.titra_point,"%5.0f") << endl;
        cout << prot.res[i_res].resID() << " " << vec2string(sum_occ,"%5.2f") << endl;
    }
    
    return 0;
}

