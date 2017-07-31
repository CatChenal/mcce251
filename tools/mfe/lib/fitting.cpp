#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "mcce.hpp"

int main(int argc, char** argv) {
    Pka pka;
    Prot prot;
    string confID;
    float shift_pK = 0.0, shift_n = 0.0;
    float threshold = 0.0;
    
    /* parsing command */
    if (argc < 2) {
        cout << "Usage:\n";
        cout << "   " << argv[0] << " confID\n";
        return -1;
    }
    
    prot.load_conflist("head3.lst");
    prot.load_occtable("fort.38");
    prot.index_on();
    vector<float> total_sum_occ(prot.titra_point.size(),0.0);

    cout << "Caution! this program prints out the occupancy of the conformer type you choose, NOT sum_crg\n" ;
    cout << "         be careful with which conformer you choose to do the titration fitting on!\n";
    
    string test_string = argv[1];
    int arg_1st_conf = 1;
    if (test_string == "-t") {
        threshold = atof(argv[2]);
        arg_1st_conf = 3;
    }
    for (int i_arg = arg_1st_conf; i_arg<argc; ++i_arg) {
        /* fitting on this conformer */
        confID = argv[i_arg];
        while (confID.size() <5) confID += ".";
        if (confID.size() > 5) while (confID.size() < 10) confID += ".";
        
        for (int i_res=0; i_res<prot.res.size(); ++i_res) {
            if (!equ_str(prot.res[i_res].resName, confID.substr(0,3))
                && !equ_str(prot.res[i_res].resID().substr(4,5), confID.substr(5,5))) continue;
            
            vector<float> sum_occ(prot.titra_point.size());
            bool match_conf = false;
            for (int i_conf=1; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
                if (equ_str(prot.res[i_res].conf[i_conf].uniqID.substr(0,confID.size()), confID)) {
                    sum_occ = sum_2vec(sum_occ, prot.res[i_res].conf[i_conf].occ_table);
                    match_conf = true;
                }
            }
            if (!match_conf) continue;
            total_sum_occ = sum_2vec(total_sum_occ, sum_occ);

            if (max(sum_occ) < threshold) continue;
            
            if (prot.param.titra_type == PH)
                cout << " ph        " << vec2string(prot.titra_point,"%5.1f") << endl;
            else
                cout << " eh        " << vec2string(prot.titra_point,"%5.0f") << endl;
            
            cout << prot.res[i_res].resID() << " " << vec2string(sum_occ,"%5.2f") << endl;
            
            pka = curve_fitting(prot.titra_point, sum_occ, shift_pK, shift_n);
            if (pka.valid) {
                if (prot.param.titra_type == EH) pka.n *= 58.;
                printf("%s %8.3f %8.3f %8.3f",prot.res[i_res].resID().c_str(), pka.pK, pka.n, 1000.*pka.chi2);
                if (pka.out_of_range) cout << " Warning! Titration may not be complete";
                cout << "\n";
            }
            else {
                printf("%s     Titration is out of range.\n",prot.res[i_res].resID().c_str());
            }
        }
    }
    if (prot.param.titra_type == PH)
        cout << " ph        " << vec2string(prot.titra_point,"%5.1f") << endl;
    else
        cout << " eh        " << vec2string(prot.titra_point,"%5.0f") << endl;
    
    cout << "Total_occ  " << vec2string(total_sum_occ,"%5.2f") << endl;
    pka = curve_fitting(prot.titra_point, total_sum_occ, shift_pK, shift_n);
    if (pka.valid) {
        if (prot.param.titra_type == EH) pka.n *= 58.;
        printf("Total_occ  %8.3f %8.3f %8.3f", pka.pK, pka.n, 1000.*pka.chi2);
        if (pka.out_of_range) cout << " Warning! Titration may not be complete";
        cout << "\n";
    }
    else {
        printf("Total_crg        Titration is out of range.\n");
    }

    cout << "Caution! this program prints out the occupancy of the conformer type you choose, NOT sum_crg\n" ;
    cout << "         be careful with which conformer you choose to do the titration fitting on!\n";
    return 0;
}


