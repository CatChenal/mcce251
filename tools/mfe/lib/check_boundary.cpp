#include "mcce.hpp"
#include "math3d.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace Math3D;

int Monte_load_pairwise(Prot* prot_p, char *dir);

int main(int argc, char **argv)
{
    Prot prot_quick, prot_default;
    int i_titra;
    float ph, eh, titra_pnt;
    float rmsd;
    int i_res, j_res, i_conf, i, i_cycle, n_cycle, n_cycle_min, n_cycle_max,n_cycle_chk;
    int i_red;
    int ic, jc;
    FILE *fp;
    time_t   timer_start, timer_end;
    
    if (argc < 3) {
        cout << argv[0] << " quick_run_dir default_run_dir" << endl;
        return -1;
    }
    
    /* quick run */
    /* load run.prm */
    string prm = argv[1];
    prm += "/run.prm";
    prot_quick.param.load_env(prm.c_str());
    prot_quick.param.scale_vdw = 0.;
    
    /* Load head3.lst */
    string conflist3 = argv[1];
    conflist3 += "/";
    conflist3 += FN_CONFLIST3;
    try {
        prot_quick.load_conflist(conflist3.c_str());
    }
    catch (string msg) {    /* if there is error in loading conflist */
        return USERERR;
    }
    if (!prot_quick.res.size()) return USERERR;
    
    prot_quick.index_on();    /* turn on extra index information */
    
    /* default run */
    /* load run.prm */
    prm = argv[2];
    prm += "/run.prm";
    prot_default.param.load_env(prm.c_str());
    prot_default.param.scale_vdw = 0.;

    /* Load head3.lst */
    conflist3 = argv[2];
    conflist3 += "/";
    conflist3 += FN_CONFLIST3;
    try {
        prot_default.load_conflist(conflist3.c_str());
    }
    catch (string msg) {    /* if there is error in loading conflist */
        return USERERR;
    }
    if (!prot_default.res.size()) return USERERR;
    
    prot_default.index_on();    /* turn on extra index information */
    
    if (prot_quick.res.size() != prot_default.res.size()) {
        cout << "Number of residues don't match between the two directories" << endl;
        return -1;
    }
    
    /* Load pairwise */
    if (Monte_load_pairwise(&prot_quick, argv[1])) return USERERR;
    if (Monte_load_pairwise(&prot_default, argv[2])) return USERERR;
    
    for (i_res = 0; i_res < prot_quick.res.size(); ++i_res) {
        int i_conf_quick, j_conf_quick, i_conf_default, j_conf_default, ic_quick, jc_quick, ic_default, jc_default;
        /* find the first charged conformer in quick run */
        for (i_conf_quick = 1; i_conf_quick < prot_quick.res[i_res].conf.size(); ++i_conf_quick) {
            if (fabs(prot_quick.res[i_res].conf[i_conf_quick].netcrg) > 0.1) break;
        }
        if (i_conf_quick >= prot_quick.res[i_res].conf.size()) continue;
        ic_quick = prot_quick.res[i_res].conf[i_conf_quick].i_sidechain_prot;
        
        /* loop over j_res */
        for (j_res = i_res+1; j_res < prot_quick.res.size(); ++j_res) {
        
            /* find the first charged conformer */
            for (j_conf_quick = 1; j_conf_quick < prot_quick.res[j_res].conf.size(); ++j_conf_quick) {
                if (fabs(prot_quick.res[j_res].conf[j_conf_quick].netcrg) > 0.1) break;
            }
            if (j_conf_quick >= prot_quick.res[j_res].conf.size()) continue;
            jc_quick = prot_quick.res[j_res].conf[j_conf_quick].i_sidechain_prot;

            /* now, find the first charged conformer in the default run */
            for (i_conf_default = 1; i_conf_default < prot_default.res[i_res].conf.size(); ++i_conf_default) {
                if (fabs(prot_default.res[i_res].conf[i_conf_default].netcrg) > 0.1) break;
            }
            if (i_conf_default >= prot_default.res[i_res].conf.size()) {
                cout << "Error! Residues don't match" << endl;
                cout << "No charged conformer found in the default run for the residue " << prot_default.res[i_res].resID() << endl;
                cout << "while there is a charged conformer in the quick run for the residue "  << prot_quick.res[i_res].resID() << " in the same slot" << endl;
                continue;
            }
            ic_default = prot_default.res[i_res].conf[i_conf_default].i_sidechain_prot;

            for (j_conf_default = 1; j_conf_default < prot_default.res[j_res].conf.size(); ++j_conf_default) {
                if (fabs(prot_default.res[j_res].conf[j_conf_default].netcrg) > 0.1) break;
            }
            if (j_conf_default >= prot_default.res[j_res].conf.size()) {
                cout << "Error! Residues don't match" << endl;
                cout << "No charged conformer found in the default run for the residue " << prot_default.res[j_res].resID() << endl;
                cout << "while there is a charged conformer in the quick run for the residue "  << prot_quick.res[j_res].resID() << " in the same slot" << endl;
                continue;
            }
            jc_default = prot_default.res[j_res].conf[j_conf_default].i_sidechain_prot;

            if (fabs(prot_quick._pairwise[ic_quick][jc_quick]) > fabs(prot_default._pairwise[ic_default][jc_default])) {
                cout << prot_quick.res[i_res].conf[i_conf_quick].uniqID << " - " << prot_quick.res[j_res].conf[j_conf_quick].uniqID << " ";
                printf("%8.3f    ", prot_quick._pairwise[ic_quick][jc_quick]);

                cout << prot_default.res[i_res].conf[i_conf_default].uniqID << " - " << prot_default.res[j_res].conf[j_conf_default].uniqID << " ";
                printf("%8.3f\n", prot_default._pairwise[ic_default][jc_default]);
            }
        }
    }
    return 0;
}

int Monte_load_pairwise(Prot* prot_p, char *dir)
{
    int   ic, jc;
    string  file_name;
    string  sbuff;
    string  stemp;
    int   i_res,i_conf,j_res,j_conf;
    
    /* declare memory */
    prot_p->_pairwise.resize(prot_p->sidechain.size());
    for (ic=0; ic<prot_p->sidechain.size(); ic++) {
        prot_p->_pairwise[ic].resize(prot_p->sidechain.size(),0.0);
    }
    
    for (ic=0; ic<prot_p->sidechain.size(); ic++) {
        if (prot_p->sidechain[ic]->uniqID.substr(3,2) == "DM") {
            continue;
        }
        
        string file_name = dir;
        file_name += "/";
        file_name += ENERGIES;
        file_name += "/" + prot_p->sidechain[ic]->uniqID + ".opp";
        load_opp(prot_p->_pairwise[ic], file_name.c_str(), *prot_p, ic);
    }
    
    return 0;
}

