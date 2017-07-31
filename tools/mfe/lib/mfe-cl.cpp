#include "mcce.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

vector < vector<double> > pairwise;
//float threshold = 0.5;

int main(int argc, char **argv) {
    Prot prot;
    ifstream fort10, conflist, fort38;
    string conflist_fname, oppfile_fname, fort38_fname;
    string uniqID, confName;
    float pairwise_threshold = 10.;
    float min_pairwise;
    int k_column = 8;
    /* read in parameters used for mfe++ */
    conflist_fname = "head3.lst";
    fort38_fname = "fort.38";
    confName = "_CL-1";
    prot.param.scale_ele = 1.;
    prot.param.scale_vdw = 1.;
    
    for (int i_arg=1; i_arg<argc; ++i_arg) {
        string flag(argv[i_arg]);
        if (flag == "-c") {
            k_column = atoi(argv[i_arg+1]);
        }
        if (flag == "-i") {
            confName = argv[i_arg+1];
        }
        if (flag == "-s") {
            prot.param.scale_vdw = atof(argv[i_arg+1]);
        }
        if (flag == "-h") {
            cout << "   Usage: mfe-cl [-c column] [-s vdw_scale] [-i confName]\n";
            return -1;
        }
    }
    
    /* open fort.38 */
    fort38.open(fort38_fname.c_str(), ios::in);
    if (!fort38) {
        cerr << "Error: can't open " << fort38_fname << " for reading. " << endl;
        return -1;
    }
    
    /* load fort.38 */
    prot.load_occtable(fort38);
    fort38.close();
    
    prot.index_on();
    pairwise.resize(prot.sidechain.size());
    for (int kc = 0; kc < prot.sidechain.size(); kc++) {
        if (prot.sidechain[kc]->confName != confName) continue;
        uniqID = prot.sidechain[kc]->uniqID;
        oppfile_fname = "energies/" + uniqID + ".opp";
        
        /* load opp file */
        ifstream  oppfile;
        oppfile.open(oppfile_fname.c_str(), ios::in);
        if (!oppfile) {
            cerr << "Error: can't open " << oppfile_fname << " for reading. " << endl;
            return -1;
        }
        
        pairwise[kc].resize(prot.sidechain.size());
        load_opp(pairwise[kc], oppfile, prot, kc);
        oppfile.close();
        
        /* open residue mfe output for writing */
        string outfilename = "confs/" + uniqID + ".out";
        FILE *outfile;
        outfile = fopen(outfilename.c_str(),"w");
        if (!outfile) {cerr << "Fail to open output file for writing\n"; exit(-1);}
        
        /* calculate mfe */
        Conf mfe_conformer;
        mfe_conformer.initialize_by_uniqID(uniqID);
        vector<float> total_mfe(prot.titra_point.size());
        
        if (prot.param.titra_type == PH) {
            fprintf(outfile, " ph        ");
            for (int i_titra=0; i_titra<prot.titra_point.size(); i_titra++) {
                fprintf(outfile, "%6.1f", prot.titra_point[i_titra]);
            }
        }
        else {
            fprintf(outfile, " eh        ");
            for (int i_titra=0; i_titra<prot.titra_point.size(); i_titra++) {
                fprintf(outfile, "%6.0f", prot.titra_point[i_titra]);
            }
        }
        fprintf(outfile, "\n");
        
        /* loop over residues */
        for (int i_res = 0; i_res < prot.res.size(); ++i_res) {
            /* define a vector to store residue mfe values */
            vector<float> mfe(prot.titra_point.size());
            /* define a vector to store total occupancy value that is counted in mfe */
            vector<float> occ_counted(prot.titra_point.size());
            
            if (!mfe_conformer.inRes(prot.res[i_res])) {    /* calculate mfe on the other residues, skip its own residue */
                if (prot.res[i_res].resName == "_CL") continue;
                if (prot.res[i_res].resName == "_CA") continue;
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
                            //if (pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot] < pairwise_threshold + min_pairwise) {
                                /* equation for mfe calculation */
                                mfe[i_titra] += pairwise[kc][prot.res[i_res].conf[i_conf].i_sidechain_prot]/1.36 * prot.res[i_res].conf[i_conf].occ_table[i_titra];
                                occ_counted[i_titra] += prot.res[i_res].conf[i_conf].occ_table[i_titra];
                            //}
                        }
                        
                        /* renormalize mfe */
                        if(occ_counted[i_titra] > 1e-5) {
                            mfe[i_titra] = mfe[i_titra] / occ_counted[i_titra];
                        }
                        else {
                            mfe[i_titra] = min_pairwise/1.36;
                        }
                    }
                    
                    total_mfe = sum_2vec(total_mfe, mfe);
                }
                
            }
            
            fprintf(outfile, "%s ", prot.res[i_res].resID().c_str());
            for (int i_titra=0; i_titra<mfe.size(); i_titra++) {
                fprintf(outfile, "%6.2f", mfe[i_titra]);
            }
            fprintf(outfile, "\n");
            
        }
        
        fprintf(outfile, "SUM        ");
        for (int i_titra=0; i_titra<prot.titra_point.size(); i_titra++) {
            fprintf(outfile, "%6.1f", total_mfe[i_titra]);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
        
        //float color = fabs(total_mfe[10]*total_mfe[10])/100.;
        //if (color > 1) color = 1;
        //if (total_mfe[10] >= 0.) {
            //printf("set_color color_%s,[1.00, %4.2f, %4.2f]\n", uniqID.substr(5,5).c_str(), 1.-color, 1.-color);
        //}
        //else {
            //printf("set_color color_%s,[%4.2f, %4.2f, 1.00]\n", uniqID.substr(5,5).c_str(), 1.-color, 1.-color);
        //}
        //printf("color color_%s, chain Y and resi %d\n", uniqID.substr(5,5).c_str(), prot.sidechain[kc]->resSeq);
        printf("Y%04d %7.2f\n", prot.sidechain[kc]->resSeq, total_mfe[k_column-1]);
    }
}
