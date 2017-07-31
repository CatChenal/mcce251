#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

vector < vector<double> > pairwise;

int main(int argc, char **argv) {
    Prot prot;
    int kc;
    ifstream conflist;
    string oppfile_fname;
    
    if (argc < 3) {
        cout << "Usage:" << endl;
        cout << "   shift_ligand residue n_ligand\n";
        return -1;
    }
    string resID = argv[1];
    int n_lig = atoi(argv[2]);
    
    /* load head3.lst and opp files */
    prot.load_conflist("head3.lst");
    prot.index_on();
    pairwise.resize(prot.sidechain.size());
    
    /* load opp files */
    prot.scale_ele = 1.;
    prot.scale_vdw = 0;
    for (kc = 0; kc < prot.sidechain.size(); kc++) {
        oppfile_fname = "energies/" + prot.sidechain[kc]->uniqID + ".opp";
        pairwise[kc].resize(prot.sidechain.size());
        
        try {
            load_opp(pairwise[kc], oppfile_fname.c_str(), prot);
        }
        catch (string& msg) {
            cout << "Assume " << prot.sidechain[kc]->uniqID << " is a dummy conformer\n";
        }
    }
    
    /* zero all pairwise within same residue and average two sides of the matrix */
    for (int ic = 0; ic < prot.sidechain.size(); ++ic) {
        for (int jc = ic; jc < prot.sidechain.size(); ++jc) {
            if (prot.sidechain[ic]->res == prot.sidechain[jc]->res) {
                pairwise[ic][jc] = pairwise[jc][ic] = 0.0;
            }
            else {
                pairwise[ic][jc] = pairwise[jc][ic] = 
                (pairwise[ic][jc] + pairwise[jc][ic]) / 2.;
            }
        }
    }
    
    int start_ires = 0;
    while (1) {
        /* find cofactor */
        int cofactor_ires, cofactor_iconf, cofactor_ic;
        for (cofactor_ires = start_ires; cofactor_ires < prot.res.size(); ++cofactor_ires) {
            if (resID.size() == 3) {
                if (prot.res[cofactor_ires].resName == resID) {
                    cout << "### Cofactor " << prot.res[cofactor_ires].resID() << endl;
                    break;
                }
            }
            else {
                if (prot.res[cofactor_ires].resName == resID.substr(0,3) &&
                prot.res[cofactor_ires].resID().substr(4,6) == resID.substr(4,6) ) {
                    cout << "### Cofactor " << prot.res[cofactor_ires].resID() << endl;
                    break;
                }
            }
        }
        if (cofactor_ires >= prot.res.size()) {
            cout << "### Done." << endl;
            break;
        }
        
        float min_pairwise = 0;
        for (int i_conf = 1; i_conf < prot.res[cofactor_ires].conf.size(); ++i_conf) {
            if (min_pairwise > get_min(pairwise[prot.res[cofactor_ires].conf[i_conf].i_sidechain_prot])) {
                min_pairwise = get_min(pairwise[prot.res[cofactor_ires].conf[i_conf].i_sidechain_prot]);
                cofactor_iconf = i_conf;
                cofactor_ic    = prot.res[cofactor_ires].conf[i_conf].i_sidechain_prot;
            }
        }
        
        /* find ligand */
        vector <int> ligand_ires(n_lig);
        vector <int> ligand_ic(n_lig);
        for (int i_lig = 0; i_lig < n_lig; ++i_lig) {
            int i_min;
            while (1) {
                string input;
                i_min = get_min_index(pairwise[cofactor_ic]);
                ligand_ires[i_lig] = prot.sidechain[i_min]->res->i_res_prot;
                ligand_ic[i_lig] = prot.sidechain[i_min]->i_sidechain_prot;
                cout << "### Ligand " << prot.sidechain[i_min]->uniqID << "? [y/n]";
                cin >> input;
                if (input == "y") break;
                else {
                    cout << "### Skip this Conformer or the whole Residue? [c/r]";
                    cin >> input;
                    if (input == "r") {
                        for (int i_conf = 1; i_conf<prot.sidechain[i_min]->res->conf.size(); ++i_conf) {
                            pairwise[cofactor_ic][prot.sidechain[i_min]->res->conf[i_conf].i_sidechain_prot] = 0.0;
                        }
                    }
                    else {
                        pairwise[cofactor_ic][i_min] = 0.0;
                    }
                }
            }
            
            /* shift all the conformers of the COFACTOR by interaction to the ligand */
            for (int i_conf = 1; i_conf < prot.res[cofactor_ires].conf.size(); ++i_conf) {
                int ic = prot.res[cofactor_ires].conf[i_conf].i_sidechain_prot;
                int jc = ligand_ic[i_lig];
                prot.res[cofactor_ires].conf[i_conf].E_extra -= pairwise[ic][jc];
            }
            
            /* reset all the conformers of the LIGAND to 0, for searching the next ligand */
            for (int i_conf = 1; i_conf < prot.res[ligand_ires[i_lig]].conf.size(); ++i_conf) {
                int ic = cofactor_ic;
                int jc = prot.res[ligand_ires[i_lig]].conf[i_conf].i_sidechain_prot;
                pairwise[ic][jc] = pairwise[jc][ic] = 0.0;
            }
        }
        
        start_ires = cofactor_ires+1;
    }
    
    /* write out new head3.lst */
    prot.write_conflist("head3.lst.new");
    
    return 0;

}

