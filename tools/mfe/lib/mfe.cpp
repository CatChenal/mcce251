#include "mcce.hpp"
#include "stdlib.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

vector < vector<double> > pairwise;
float threshold = 0.5;

int main(int argc, char **argv) {
    Prot prot;
    ifstream fort10, conflist, oppfile, fort38;
    string conflist_fname, oppfile_fname, fort38_fname;
    string uniqID;
    float pairwise_threshold = 10.;
    float min_pairwise;
    const char* USAGE="\
                Usage:\n\
                      mfe++ conformer_name [-t threshold_for_print] [-x threshold_for_exclusion]\n\
                      You can also run mfe++ without parameter, and using fort.10 instead:\n\
                      1. put the path to the head3.lst in the first line of fort.10;\n\
                      2. put the path to the opp file in the second line of fort.10;\n\
                      3. put the path to the occupany file (fort.38) in the third line of fort.10;\n\
                      4. put the conformer name in the forth line of fort.10;\n\
                      5. run mfe++\n\
                      The output is a tab delimited file (Cat's version, Jan 2014)\n";

    /* read in parameters used for mfe++ */
    if (argc == 1) {
        /* read fort.10 */
        fort10.open("fort.10", ios::in);
        if (!fort10) {
            cout << USAGE << endl;
            return -1;
        }
        if ( !getline(fort10, conflist_fname) || !getline(fort10, oppfile_fname) || !getline(fort10, fort38_fname) || !getline(fort10, uniqID) ) {
            cerr << "Error: fort.10 not in correct format" << endl;
            cout << USAGE << endl;
            return -1;
        }
        fort10.close();
    }
    else {
        conflist_fname = "head3.lst";
        uniqID = argv[1];
        oppfile_fname = "energies/" + uniqID + ".opp";
        fort38_fname = "fort.38";
        for (int i_arg=2; i_arg<argc-1; ++i_arg) {
            string flag(argv[i_arg]);
            if (flag == "-t") {
                threshold = atof(argv[i_arg+1]);
            }
            if (flag == "-x") { /* exclude */
                pairwise_threshold=atof(argv[i_arg+1]);
            }
        }
    }
    /* open fort.38 */
    fort38.open(fort38_fname.c_str(), ios::in);
    if (!fort38) {
        cerr << "Error: can't open " << fort38_fname << " for reading. " << endl;
        return -1;
    }
    /* open opp file */
    oppfile.open(oppfile_fname.c_str(), ios::in);
    if (!oppfile) {
        cerr << "Error: can't open " << oppfile_fname << " for reading. " << endl;
        return -1;
    }
    /* load fort.38 */
    prot.load_occtable(fort38);
    fort38.close();
    /* load opp file */
    prot.index_on();
    pairwise.resize(prot.sidechain.size());
    int kc;
    for (kc = 0; kc < prot.sidechain.size(); kc++) {
        if (prot.sidechain[kc]->uniqID == uniqID) {
            break;
        }
    }
    if (kc >= prot.sidechain.size()) {
        cerr << "Fail to find conformer " << uniqID << " in " << fort38_fname << ' ' << prot.res.size() << endl;
        return -1;
    }
    prot.param.scale_ele = 1.;
    prot.param.scale_vdw = 1.;

    pairwise[kc].resize(prot.sidechain.size());
    load_opp(pairwise[kc], oppfile, prot, kc);
    oppfile.close();
    /* open residue mfe output for writing */
    FILE *outfile;
    string filename = uniqID;
    filename += ".mfe";
    outfile = fopen(filename.c_str(),"w");
    if (!outfile) {cerr << "Fail to open " << filename << " for writing\n"; exit(-1);}

    /* calculate mfe */
    Conf mfe_conformer;
    mfe_conformer.initialize_by_uniqID(uniqID);
    vector<float> total_mfe(prot.titra_point.size());

    if (prot.param.titra_type == PH) {
        printf("ph        \t");
        fprintf(outfile, "ph        \t");
        for (int i_titra=0; i_titra<prot.titra_point.size(); i_titra++) {
            printf("%6.1f\t", prot.titra_point[i_titra]);
            fprintf(outfile, "%6.1f\t", prot.titra_point[i_titra]);
        }
    }
    else {
        printf("eh        \t");
        fprintf(outfile, "eh        \t");
        for (int i_titra=0; i_titra<prot.titra_point.size(); i_titra++) {
            printf("%6.0f\t", prot.titra_point[i_titra]);
            fprintf(outfile, "%6.0f\t", prot.titra_point[i_titra]);
        }
    }
    printf("\n");
    fprintf(outfile, "\n");

    /* loop over residues */
    for (int i_res = 0; i_res < prot.res.size(); ++i_res) {
        /* define a vector to store residue mfe values */
        vector<float> mfe(prot.titra_point.size());
        /* define a vector to store total occupancy value that is counted in mfe */
        vector<float> occ_counted(prot.titra_point.size());

        if (!mfe_conformer.inRes(prot.res[i_res])) {    /* calculate mfe on the other residues, skip its own residue */
                if (prot.res[i_res].conf.size() > 1) {  /* calculate mfe on residues that have sidechain conformers */

                /* loop over all titration points */
                for (int i_titra = 0; i_titra<prot.res[i_res].conf[1].occ_table.size(); i_titra++) {
                    int i_conf;
                    /* initialize min_pairwise between mfe_conformer and i_res */
                    for (i_conf = 1; i_conf < prot.res[i_res].conf.size(); ++i_conf) {
                        if (prot.res[i_res].conf[i_conf].occ_table[i_titra] < 1e-5) continue;
                        min_pairwise = pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot];
                        break;
                    }
                    if (i_conf >= prot.res[i_res].conf.size()) { /* no conformer occupied */
                        min_pairwise = 0.;
                    }
                    /* search for the minimum value */
                    for (i_conf = 1; i_conf < prot.res[i_res].conf.size(); ++i_conf) {
                        if (prot.res[i_res].conf[i_conf].occ_table[i_titra] < 1e-5) continue;
                        if (pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot]<min_pairwise)
                            min_pairwise=pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot];
                    }
                    /* calculate mfe */
                    for (i_conf = 1; i_conf < prot.res[i_res].conf.size(); ++i_conf) {
                        /* equation for mfe calculation */
                        mfe[i_titra] += pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot]/PH2KCAL * prot.res[i_res].conf[i_conf].occ_table[i_titra];
                        occ_counted[i_titra] += prot.res[i_res].conf[i_conf].occ_table[i_titra];
                    }
                    /* renormalize mfe */
                    if(occ_counted[i_titra] > 1e-5) {
                        mfe[i_titra] = mfe[i_titra] / occ_counted[i_titra];
                    }
                    else {
                        mfe[i_titra] = min_pairwise/PH2KCAL;
                    }
                }
                total_mfe = sum_2vec(total_mfe, mfe);
            }
        }
/*175 */
        fprintf(outfile, "%s\t", prot.res[i_res].resID().c_str());
        for (int i_titra=0; i_titra<mfe.size(); i_titra++) {
            fprintf(outfile, "%5.1f\t", mfe[i_titra]);
        }
        fprintf(outfile, "\n");

        if (max_abs(mfe) > threshold) {
            printf("%s\t", prot.res[i_res].resID().c_str());
            for (int i_titra=0; i_titra<mfe.size(); i_titra++) {
                printf("%5.1f\t", mfe[i_titra]);
            }
            printf("\n");
        }
    }
    printf("SUM      \t");
    fprintf(outfile, "SUM      \t");
    for (int i_titra=0; i_titra<prot.titra_point.size(); i_titra++) {
        printf("%6.1f\t", total_mfe[i_titra]);
        fprintf(outfile, "%6.1f\t", total_mfe[i_titra]);
    }
    printf("\n");
    fprintf(outfile, "\n");
}
