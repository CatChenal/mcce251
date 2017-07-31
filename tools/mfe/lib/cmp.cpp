#include <string>
#include <vector>
#include "mcce.hpp"

using namespace std;

int main(int argc, char** argv) {
    Prot prot1, prot2;
    float threshold = 0.05;
    
    if (argc < 3) {
        cout << argv[0] << " dir1 dir2\n";
        return -1;
    }
    if (argc >= 3) {
        for (int iarg=3; iarg<argc-1; ++iarg) {
            string arg_str(argv[iarg]);
            if (arg_str == "-t") { threshold = atof(argv[iarg+1]); }
        }
    }
    
    string dir1(argv[1]), dir2(argv[2]);
    string filename;
    
    filename = dir1 + "/head3.lst";
    prot1.load_conflist(filename.c_str());
    filename = dir1 + "/fort.38";
    prot1.load_occtable(filename.c_str());
    prot1.make_sum_crg();
    
    filename = dir2 + "/head3.lst";
    prot2.load_conflist(filename.c_str());    
    filename = dir2 + "/fort.38";
    prot2.load_occtable(filename.c_str());
    prot2.make_sum_crg();
    
    int n_titra = prot1.titra_point.size() < prot2.titra_point.size() ?
    prot1.titra_point.size():prot2.titra_point.size();
    
    for (int i_res =0 ; i_res<prot1.res.size(); ++i_res) {
        Res* res_p = prot2.find_res(prot1.res[i_res]);
        if (!res_p) continue;
        
        vector <float> diff(n_titra);
        for (int i_titra=0; i_titra<n_titra; ++i_titra) {
            diff[i_titra] = prot1.res[i_res].sum_H[i_titra] - prot2.res[i_res].sum_H[i_titra];
        }
        if (max_abs(diff) > threshold) {
            cout << prot1.res[i_res];
            print_vec(diff, "%5.2f");
        }
    }
}

