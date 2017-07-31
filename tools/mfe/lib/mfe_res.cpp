#include "mcce.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

extern vector < vector<double> > pairwise;
extern vector<Res*>            flip_res;   /* a list of residue to flip in each MC step */
extern long        dummy_num;      /* seed for random number generator */
extern double      beta;           /* = 1/kT */
float scale_ele = 1.;
float scale_vdw = 0.25;
float ph = 7.0;
float eh = 0.0;
extern void Monte_MC(Prot&);
extern Prot Monte_reduce(const Prot& prot);

int main(int argc, char **argv) {
    Prot prot;
    vector<string> search_list;
    
    if (argc < 2) {
        cout << "Usage:" << endl;
        cout << "   " << argv[0] << " resID\n";
        return -1;
    }
    
    for (int i_arg = 1; i_arg<argc; ++i_arg) {
        search_list.push_back(argv[i_arg]);
        while (search_list.back().size() < 10) {
            search_list.back().push_back('.');
        }
    }
    
    /* load run.prm */
    prot.param.load_env("run.prm");
    
    /* load head3.lst and opp files */
    prot.load_conflist("head3.lst");
    prot.load_occtable("fort.38");
    prot.index_on();
    pairwise.resize(prot.sidechain.size());
    for (int ic=0; ic<prot.sidechain.size(); ic++) {
        pairwise[ic].resize(prot.sidechain.size(),0.0);
    }

    /* scale energy */
    prot.param.scale_ele     = 1.0;
    prot.param.scale_dsolv   = 1.0;
    prot.param.scale_tor     = 0.5;
    prot.param.scale_vdw0    = 0.25;
    prot.param.scale_vdw1    = 0.25;
    prot.param.scale_vdw     = 0.25;
    for (int ic=0; ic<prot.sidechain.size(); ic++) {
        prot.sidechain[ic]->E_vdw0    *= prot.param.scale_vdw0;
        prot.sidechain[ic]->E_vdw1    *= prot.param.scale_vdw1;
        prot.sidechain[ic]->E_epol    *= prot.param.scale_ele;
        prot.sidechain[ic]->E_tors    *= prot.param.scale_tor;
        prot.sidechain[ic]->E_dsolv   *= prot.param.scale_dsolv;
    }
    
    /* turn off all conformers */
    for (int ic=0; ic<prot.sidechain.size(); ic++) {
        prot.sidechain[ic]->toggle = 'f';
    }
    
    /* turn on those residue that match search list and load opp for them */
    for (int i_sch = 0; i_sch< search_list.size(); ++i_sch) {
        /* find the matching residue */
        int k_res;
        for (k_res=0; k_res<prot.res.size(); ++k_res) {
            if (equ_str(search_list[i_sch],prot.res[k_res].resID())) break;
        }
        if (k_res >= prot.res.size()) {
            cerr << "Cannot find residue " << search_list[i_sch] << endl;
            return -1;
        }
        for (int k_conf=1; k_conf<prot.res[k_res].conf.size(); ++k_conf) {
            prot.res[k_res].conf[k_conf].toggle = 't';
            
            if (prot.res[k_res].conf[k_conf].uniqID.substr(3,2) == "DM") {
                continue;
            }
            
            string file_name = ENERGIES;
            file_name += "/" + prot.res[k_res].conf[k_conf].uniqID + ".opp";
            int ic = prot.res[k_res].conf[k_conf].i_sidechain_prot;
            load_opp(pairwise[ic], file_name.c_str(), prot);
        }
    }
    
    /* use n_flip in run.prm to initialize flip_res array */
    prot.param.monte_flips = 3;
    flip_res.resize(prot.param.monte_flips);
    dummy_num = time(NULL);
    prot.param.monte_temp = 273.15;
    beta = KCAL2KT/(prot.param.monte_temp/ROOMT);

    for (int i_titra = 0; i_titra<prot.titra_point.size(); i_titra++) {
        
        /* calc. pH and Eh dependent self-energy of each conformer */
        for (int ic=0; ic<prot.sidechain.size(); ic++) {
            prot.sidechain[ic]->E_ph =
            prot.param.monte_temp/ROOMT * prot.sidechain[ic]->H * (ph-prot.sidechain[ic]->pKa) * PH2KCAL;
            
            prot.sidechain[ic]->E_eh =
            prot.param.monte_temp/ROOMT * prot.sidechain[ic]->e * (eh-prot.sidechain[ic]->Em) * PH2KCAL/58.0;
            
            prot.sidechain[ic]->E_self
            = prot.sidechain[ic]->E_vdw0
            + prot.sidechain[ic]->E_vdw1
            + prot.sidechain[ic]->E_epol
            + prot.sidechain[ic]->E_tors
            + prot.sidechain[ic]->E_dsolv
            + prot.sidechain[ic]->E_extra
            + prot.sidechain[ic]->E_ph
            + prot.sidechain[ic]->E_eh;
        }
        
        Prot prot_w = Monte_reduce(prot);
        prot_w.index_off();
        prot_w.index_on();

        prot_w.weighted_res_list.clear();
        for (int i_res=0;i_res<prot_w.res.size();++i_res) {
            prot_w.weighted_res_list.push_back(i_res);
        }
        
        Monte_MC(prot_w);
        for (int ic=0;ic<prot_w.sidechain.size();ic++) {
            prot_w.sidechain[ic]->occ_table.push_back(prot_w.sidechain[ic]->occ);
        }
        
        if (prot_w.param.titra_type == PH) {
            printf(" ph           %s\n",vec2string(prot_w.titra_point,"%5.1f").c_str());
        }
        else {
            printf(" eh           %s\n",vec2string(prot_w.titra_point,"%5.0f").c_str());
        }
        for (int ic=0; ic<prot_w.sidechain.size(); ic++) {
            printf("%s %s\n",prot_w.sidechain[ic]->uniqID.c_str(), vec2string(prot_w.sidechain[ic]->occ_table, "%5.3f").c_str());
        }
    }
    
}

/*
vector< vector< vector<float> > > conf_to_res_mfe(prot.res[k_res].conf.size());
vector< vector<float> > conf_to_prot_mfe(prot.res[k_res].conf.size());
for (int k_conf = 1; k_conf<prot.res[k_res].conf.size(); ++k_conf) {
    conf_to_res_mfe[k_conf].resize(prot.res.size());
    conf_to_prot_mfe[k_conf].resize(prot.res.size());
}

calc. MFE 
for (int k_conf = 1; k_conf<prot.res[k_res].conf.size(); ++k_conf) {
    int kc = prot.res[k_res].conf[k_conf].i_sidechain_prot;
    
    for (int i_res = 0; i_res < prot.res.size(); ++i_res) {
        conf_to_res_mfe[k_conf][i_res].resize(prot.titra_point.size(), 0.0);
        if (i_res != k_res) continue;
        
        for (int i_conf = 1; i_conf < prot.res[i_res].conf.size(); ++i_conf) {
            if (prot.param.titra_type == PH) ph = prot.titra_point[i_titra];
            else eh = prot.titra_point[i_titra];
            
            conf_to_res_mfe[k_conf][i_res][i_titra]
            += pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot]/1.36
            * prot.res[i_res].conf[i_conf].occ_table[i_titra];
            
        }
    }
    conf_to_prot_mfe[k_conf] = sum_2vec(conf_to_prot_mfe[k_conf], conf_to_res_mfe[k_conf][i_res]);
}
    }
    */
/*
vector<string>::iterator pos;
pos = find(conf_type.begin(),conf_type.end(),confName);
*/
