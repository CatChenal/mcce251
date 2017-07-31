#include "mcce.hpp"
#include <sstream>

int load_opp(vector <double>& pairwise_kconf, ifstream& oppfile, const Prot& prot, const int kconf)
{
    string line;
    int ic, ic_start;

    pairwise_kconf.clear();
    pairwise_kconf.resize(prot.sidechain.size(),0.0);

    for (ic=0; ic<prot.sidechain.size(); ic++)
        prot.sidechain[ic]->on = 0;

    ic = 0;
    while (getline(oppfile, line)) {
        line = rm_comment(line);
        if (line.size() < 20) break;

        stringstream line_ss(line);
        int serial;
        string uniqID;
        double ele_pair, vdw_pair;

        line_ss >> serial;
        line_ss >> uniqID;
        line_ss >> ele_pair;
        line_ss >> vdw_pair;

        ic_start = ic;
        while (uniqID != prot.sidechain[ic]->uniqID) {
            ic++;
            if (ic >= prot.sidechain.size()) ic = 0;
            if (ic == ic_start) break;
        }
        if (uniqID == prot.sidechain[ic]->uniqID) {
            /* rescale based on scaling parameter */
            if (vdw_pair < 900. || fabs(prot.param.scale_vdw) < 1e-4) vdw_pair *= prot.param.scale_vdw;
            if (fabs(ele_pair) < 900. || fabs(prot.param.scale_ele) < 1e-4) ele_pair *= prot.param.scale_ele;

            pairwise_kconf[ic] = ele_pair + vdw_pair;
            prot.sidechain[ic]->on = 1;
        }
    }

    /* check if all interactions are loaded */
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        if (prot.sidechain[ic]->on) continue;
        if (prot.sidechain[ic]->uniqID.substr(3,2) == "DM") continue;
        cerr << "Error: interaction to "<< prot.sidechain[ic]->uniqID <<" is not loaded for "
        << prot.sidechain[kconf]->uniqID << "\n";
    }
    return 0;
}

int load_opp(vector <double>& pairwise_kconf, const char *file_name, const Prot& prot, const int kconf)
{
    ifstream oppfile;
    string msg;

    oppfile.open(file_name, ios::in);
    if (!oppfile) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }

    load_opp(pairwise_kconf, oppfile, prot, kconf);
    oppfile.close();
    return 0;
}

int load_opp(double* pairwise_kconf, ifstream& oppfile, const Prot& prot)
{
    string line;
    int ic, ic_start;

    for (ic=0; ic<prot.sidechain.size(); ic++)
        prot.sidechain[ic]->on = 0;

    ic = 0;
    while (getline(oppfile, line)) {
        line = rm_comment(line);
        if (line.size() < 20) break;

        stringstream line_ss(line);
        int serial;
        string uniqID;
        double ele_pair, vdw_pair;

        line_ss >> serial;
        line_ss >> uniqID;
        line_ss >> ele_pair;
        line_ss >> vdw_pair;

        ic_start = ic;
        while (uniqID != prot.sidechain[ic]->uniqID) {
            ic++;
            if (ic >= prot.sidechain.size()) ic = 0;
            if (ic == ic_start) break;
        }
        if (uniqID == prot.sidechain[ic]->uniqID) {
            /* rescale based on scaling parameter */
            if (vdw_pair < 900. || prot.param.scale_vdw < 1e-4) vdw_pair *= prot.param.scale_vdw;
            if (fabs(ele_pair) < 900. || prot.param.scale_ele < 1e-4) ele_pair *= prot.param.scale_ele;

            pairwise_kconf[ic] = ele_pair + vdw_pair;
            prot.sidechain[ic]->on = 1;
        }
    }

    /* check if all interactions are loaded */
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        if (prot.sidechain[ic]->on) continue;
        if (prot.sidechain[ic]->uniqID.substr(3,2) == "DM") continue;
        cerr << "Error: interaction to "<< prot.sidechain[ic]->uniqID <<" is not loaded\n";
    }
    return 0;
}

int load_opp(double* pairwise_kconf, const char *file_name, const Prot& prot)
{
    ifstream oppfile;
    string msg;

    oppfile.open(file_name, ios::in);
    if (!oppfile) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }

    load_opp(pairwise_kconf, oppfile, prot);
    oppfile.close();
    return 0;
}

