#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

int main() {
    Prot prot;
    
    prot.load_conflist("head3.lst");
    prot.index_on();

    for (int i_res=0; i_res<prot.res.size(); ++i_res) {
        for (int i_conf=1; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
            if (prot.res[i_res].conf[i_conf].uniqID.substr(0,10) != "ARG+1A0482") continue;
            if (prot.res[i_res].conf[i_conf].uniqID != "ARG+1A0482_035" &&
                prot.res[i_res].conf[i_conf].uniqID != "ARG+1A0482_036" &&
                prot.res[i_res].conf[i_conf].uniqID != "ARG+1A0482_037" &&
                prot.res[i_res].conf[i_conf].uniqID != "ARG+1A0482_038") {
                prot.res[i_res].conf[i_conf].toggle = 't';
            }
        }
    }
    
    prot.write_conflist("head3.lst");
    return 0;
}

