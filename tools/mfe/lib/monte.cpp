/* This is Monte Carlo part of MCCE */

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "mcce.hpp"

using namespace std;

int  Monte_load_pairwise(Prot&);
void Monte_get_biglist(Prot&);
void Monte_check_toggle(Prot&);
void Monte_set_toggle(Prot&);
Prot Monte_reduce(const Prot& prot);
void Monte_MC(Prot&);
int  Monte_out(Prot& prot);
void Monte_zero_counters(Prot&);
void Monte_update_TS(Prot& prot);
void Monte_group_conf_type(Prot& prot);
void Monte_update_energy(Prot& prot);

vector< vector<double> > pairwise;  /* pairwise interaction between conf-conf */
vector<Res*>            flip_res;   /* a list of residue to flip in each MC step */
long        dummy_num;      /* seed for random number generator */
double      beta;           /* = 1/kT */

int Monte()
{
    Prot prot, prot_w, prot_red;
    int i_titra;
    float ph, eh; // titra_pnt;
    float rmsd;
    int i_res, j_res, i_conf, i, i_cycle, n_cycle, n_cycle_min, n_cycle_max,n_cycle_chk;
    int i_red;
    int ic, jc;
    FILE *fp;
    time_t   timer_start, timer_end;
    
    remove(MC_LOG);
    
    /* Record starting time */
    timer_start = time(NULL);

    /* load run.prm */
    prot.param.load_env("run.prm");
    if (prot.param.extra.size()) {
        prot.param.load_extra(prot.param.extra.c_str());
    }

    printf("   Scaling factors:\n");
    printf("   VDW0  = %.3f:\n", prot.param.scale_vdw0);
    printf("   VDW1  = %.3f:\n", prot.param.scale_vdw1);
    printf("   VDW   = %.3f:\n", prot.param.scale_vdw);
    printf("   TORS  = %.3f:\n", prot.param.scale_tor);
    printf("   ELE   = %.3f:\n", prot.param.scale_ele);
    printf("   DSOLV = %.3f:\n\n", prot.param.scale_dsolv);
    
    /* Load head3.lst */
    try {
        prot.load_conflist(FN_CONFLIST3);
    }
    catch (string msg) {    /* if there is error in loading conflist */
        return USERERR;
    }
    if (!prot.res.size()) return USERERR;
    
    prot.index_on();    /* turn on extra index information */
    
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        /* switch toggle */
        if (prot.sidechain[ic]->toggle == 't') prot.sidechain[ic]->toggle = 'f';
        else if (prot.sidechain[ic]->toggle == 'f') prot.sidechain[ic]->toggle = 't';
        
        /* scale energy */
        prot.sidechain[ic]->E_vdw0    *= prot.param.scale_vdw0;
        prot.sidechain[ic]->E_vdw1    *= prot.param.scale_vdw1;
        prot.sidechain[ic]->E_epol    *= prot.param.scale_ele;
        prot.sidechain[ic]->E_tors    *= prot.param.scale_tor;
        prot.sidechain[ic]->E_dsolv   *= prot.param.scale_dsolv;

        /* Backup i_sidechain_prot, which is needed for extracting pairwise during reduced run */
        prot.sidechain[ic]->i_sidechain_prot_orig
        = prot.sidechain[ic]->i_sidechain_prot;
    }

    /* use n_flip in run.prm to initialize flip_res array */
    flip_res.resize(prot.param.monte_flips);
    
    /* Load pairwise */
    if (Monte_load_pairwise(prot)) return USERERR;
    
    Monte_check_toggle(prot);
    
    /* initialize E_base */
    prot.E_base = 0;
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        prot.sidechain[ic]->E_base = 0;
    }
    
    if (prot.param.monte_seed < 0) dummy_num = time(NULL);
    else dummy_num = prot.param.monte_seed;
    timer_end = time(NULL);
    
    {
        fp = fopen(MC_LOG,"w"); 
        fprintf(fp,"random number seed = %ld\n",dummy_num);
        fprintf(fp,"Monte Carlo set up time: %ld seconds.\n", timer_end-timer_start); fflush(stdout);
        fprintf(fp,"Starting number of residue = %d, number of conformer = %d\n",(int)prot.res.size(),(int)prot.sidechain.size());
        fprintf(fp,"   Scaling factors:\n");
        fprintf(fp,"   VDW0  = %.3f:\n", prot.param.scale_vdw0);
        fprintf(fp,"   VDW1  = %.3f:\n", prot.param.scale_vdw1);
        fprintf(fp,"   VDW   = %.3f:\n", prot.param.scale_vdw);
        fprintf(fp,"   TORS  = %.3f:\n", prot.param.scale_tor);
        fprintf(fp,"   ELE   = %.3f:\n", prot.param.scale_ele);
        fprintf(fp,"   DSOLV = %.3f:\n", prot.param.scale_dsolv);
        fclose(fp);
    }
    
    remove(MC_DETAIL);

    for (i=0;i<prot.param.monte_nstart;i++) ran2(&dummy_num);
    
    /* pH/Eh titration */
    prot.titra_point.clear();
    for (i_titra=0;i_titra<prot.param.titra_steps;i_titra++) {
        timer_start = time(NULL);
        
        ph = prot.param.titra_ph0;
        eh = prot.param.titra_eh0;
        if (prot.param.titra_type == PH) {
            ph += ((float)i_titra)*prot.param.titra_phd;
            prot.titra_point.push_back(ph);
        }
        else {
            eh += ((float)i_titra)*prot.param.titra_ehd;
            prot.titra_point.push_back(eh);
        }        
        fp = fopen(MC_LOG,"a"); fprintf(fp,"\npH = %6.2f  Eh = %6.2f\n", ph, eh); fclose(fp);
        
        /* calc. pH and Eh dependent self-energy of each conformer
        TS is set to 0 at the beginning */
        for (ic=0; ic<prot.sidechain.size(); ic++) {
            prot.sidechain[ic]->E_ph =
            prot.param.monte_temp/ROOMT * prot.sidechain[ic]->H * (ph-prot.sidechain[ic]->pKa) * PH2KCAL;
            prot.sidechain[ic]->E_eh =
            prot.param.monte_temp/ROOMT * prot.sidechain[ic]->e * (eh-prot.sidechain[ic]->Em) * MV2KCAL;
            prot.sidechain[ic]->TS   =  0.0;
            
            prot.sidechain[ic]->E_self
            = prot.sidechain[ic]->E_base
            + prot.sidechain[ic]->E_vdw0
            + prot.sidechain[ic]->E_vdw1
            + prot.sidechain[ic]->E_epol
            + prot.sidechain[ic]->E_tors
            + prot.sidechain[ic]->E_dsolv
            + prot.sidechain[ic]->E_extra
            + prot.sidechain[ic]->E_ph
            + prot.sidechain[ic]->E_eh
            + prot.sidechain[ic]->TS;       /* entropy correction */
        }
        
        prot_w = Monte_reduce(prot);
        prot_w.index_off();
        prot_w.index_on();
        for (ic=0;ic<prot_w.sidechain.size();ic++) {
            prot_w.sidechain[ic]->occ_table.clear();
        }
        
        i_red =0;
        while(1) {
            /* one reduce run */
            
            if (prot.param.monte_n_red > 0) {
                /* reaches maximum reducing steps */
                if (i_red >= prot.param.monte_n_red) break;
            }
            if (!prot_w.sidechain.size()) break;    /* no conformer to sample */
            
            /* calculate deviation of occupancies from finished MC to decide if it is converged */
            if (i_red >= 2) {
                float max_dev = 0.;
                rmsd = 0.;
                for (ic=0;ic<prot_w.sidechain.size();ic++) {
                    float dev_occ = stdev(prot_w.sidechain[ic]->occ_table);
                    if (dev_occ > max_dev) max_dev = dev_occ;
                    rmsd += dev_occ*dev_occ;
                }
                rmsd = sqrt(rmsd/(float)prot_w.sidechain.size());
                /* information up to the last reduce run */
                fp = fopen(MC_LOG,"a"); fprintf(fp,"Reduce = %3d, rmsd = %.5f, max_dev = %.4f\n",i_red,rmsd,max_dev); fclose(fp);
                if (max_dev < prot.param.monte_converge) break;
            }
            
            if (i_red >=2) Monte_set_toggle(prot_w); /* make sure first two reduce runs are full runs */
            prot_red = Monte_reduce(prot_w);
            prot_red.index_off();
            prot_red.index_on();
            for (ic=0;ic<prot_red.sidechain.size();ic++) {
                prot_red.sidechain[ic]->occ_table.clear();
            }

            if (!prot_red.sidechain.size()) break;
            Monte_get_biglist(prot_red);
            
            fp = fopen(MC_LOG,"a"); fprintf(fp,"\nReduce = %3d, number of active residue: %4d, active conformer %5d\n",i_red+1,(int)prot_red.res.size(),(int)prot_red.sidechain.size()); fclose(fp);
            
            /* Initialize */
            for (i_res=0;i_res<prot_red.res.size();i_res++) {
                i_conf = (int) (ran2(&dummy_num) * (prot_red.res[i_res].conf.size() - 1.)) + 1;
                prot_red.res[i_res].conf_w = &prot_red.res[i_res].conf[i_conf];
            }
            Monte_zero_counters(prot_red);
            
            prot_red.E_state = prot_red.E_base;
            //printf("%f\n",prot_red.E_state);
            for (i_res=0;i_res<prot_red.res.size();i_res++) {
                prot_red.E_state += prot_red.res[i_res].conf_w->E_self;
                ic = prot_red.res[i_res].conf_w->i_sidechain_prot_orig;
                for (j_res=i_res+1;j_res<prot_red.res.size();j_res++) {
                    jc = prot_red.res[j_res].conf_w->i_sidechain_prot_orig;
                    prot_red.E_state += pairwise[ic][jc];
                }
            }
            prot_red.E_min = prot_red.E_state;

            /* Equilibrating */
            n_cycle = (prot_red.sidechain.size()*prot.param.monte_neq - 1)/prot.param.monte_niter_cycle + 1;
            beta = KCAL2KT/(prot.param.monte_temp/ROOMT);
            fp = fopen(MC_LOG,"a");
            fprintf(fp,"Equilibrating:   %5d cycle(s) at temperature = %10.2fK\n",n_cycle,prot.param.monte_temp);
            fclose(fp);
            
            /* initialize weighted residue list */
            prot_red.weighted_res_list.clear();
            for (i_res=0;i_res<prot_red.res.size();++i_res) {
                prot_red.weighted_res_list.push_back(i_res);
            }
            //cout << vec2string(prot_red.weighted_res_list, "%3d") << endl;
            for (i_cycle=0;i_cycle<n_cycle;i_cycle++) {
                Monte_MC(prot_red);
                for (ic=0;ic<prot_red.sidechain.size();ic++) {
                    prot_red.sidechain[ic]->occ = (float) prot_red.sidechain[ic]->counter_accept / (float) ((i_cycle+1) * prot.param.monte_niter_cycle);
                }
                if (i_red == 0) {
                    Monte_update_energy(prot_red); /* update TS correction in the first reduce run, later it is inherited from prot_w */
                    /* update prot_red.E_state based on res.conf_w */
                    prot_red.E_state = prot_red.E_base;
                    for (int i_res=0;i_res<prot_red.res.size();i_res++) {
                        prot_red.E_state += prot_red.res[i_res].conf_w->E_self;
                        int ic = prot_red.res[i_res].conf_w->i_sidechain_prot_orig;
                        for (int j_res=i_res+1;j_res<prot_red.res.size();j_res++) {
                            int jc = prot_red.res[j_res].conf_w->i_sidechain_prot_orig;
                            prot_red.E_state += pairwise[ic][jc];
                        }
                    }
                }
            }
            
            /* Collect Data */
            n_cycle_min = (prot_red.sidechain.size()*prot.param.monte_niter_min - 1)/prot.param.monte_niter_cycle + 1;
            if (prot.param.monte_niter_max > 0)
                n_cycle_max = (prot_red.sidechain.size()*prot.param.monte_niter_max - 1)/prot.param.monte_niter_cycle + 1;
            else n_cycle_max = -1;
            n_cycle_chk = (prot_red.sidechain.size()*prot.param.monte_niter_chk - 1)/prot.param.monte_niter_cycle + 1;
            beta = KCAL2KT/(prot.param.monte_temp/ROOMT);
            fp = fopen(MC_LOG,"a");
            fprintf(fp,"Collecting data:                at temperature = %10.2fK\n",prot.param.monte_temp);
            fprintf(fp,"   Minimum/maximum number of cycles: %5d/%5d\n",n_cycle_min,n_cycle_max);
            fclose(fp);
            
            i_cycle = 0;
            
            while (1) {
                /* one cycle */
                double E_chk = prot_red.E_base;
                for (i_res=0;i_res<prot_red.res.size();i_res++) {
                    E_chk += prot_red.res[i_res].conf_w->E_self;
                    int ic = prot_red.res[i_res].conf_w->i_sidechain_prot_orig;
                    for (j_res=i_res+1;j_res<prot_red.res.size();j_res++) {
                        int jc = prot_red.res[j_res].conf_w->i_sidechain_prot_orig;
                        E_chk += pairwise[ic][jc];
                    }
                }
                //printf("i_cycle=%3d, E_state=%10.4f,E_chk=%10.4f,diff=%.3e\n",i_cycle,prot_red.E_state,E_chk,prot_red.E_state-E_chk);
                
                
                if (i_cycle > n_cycle_min) { /* check the convergence */
                    if (n_cycle_max > 0) {
                        if (i_cycle > n_cycle_max) break;
                    }
                    
                    if (i_cycle % n_cycle_chk == 0) {
                        float max_dev = 0, occ_dev;
                        /*
                        moved to after Monte_MC() 03/02/2007
                        for (ic=0;ic<prot_red.sidechain.size();ic++) {
                            prot_red.sidechain[ic]->occ_old = prot_red.sidechain[ic]->occ;
                            prot_red.sidechain[ic]->occ = (float) prot_red.sidechain[ic]->counter_accept / (float) (i_cycle * prot.param.monte_niter_cycle);
                        }
                        */
                        
                        //rmsd = 0.;
                        for (ic=0;ic<prot_red.sidechain.size();ic++) {
                            occ_dev = stdev(prot_red.sidechain[ic]->occ_table);
                            if (occ_dev > max_dev) max_dev = occ_dev;
                        }
                        if (max_dev < prot.param.monte_converge) break;
                    }
                }
                
                /* update weighted residue list */
                prot_red.weighted_res_list.clear();
                for (i_res=0;i_res<prot_red.res.size();++i_res) {
                    vector<float> tmp_devs;
                    if (i_cycle > 1) {
                        for (i_conf = 1; i_conf<prot_red.res[i_res].conf.size(); ++i_conf) {
                            //cout << vec2string(prot_red.res[i_res].conf[i_conf].occ_table, "%5.2f") << endl;
                            tmp_devs.push_back(stdev(prot_red.res[i_res].conf[i_conf].occ_table));
                        }
                    }
                    int weight = 1;
                    //int weight = 1 + (int)(1000. * max(tmp_devs));
                    //cout << "ires= " << i_res << " weight= " << weight << endl;
                    //cout << "dev vector " << vec2string(tmp_devs, "%5.2f") << " weight " << weight << endl;
                    prot_red.weighted_res_list.resize(prot_red.weighted_res_list.size()+weight,i_res);
                }
                
                //cout << vec2string(prot_red.weighted_res_list, "%3d") << endl;
                
                /* Set zero */
                Monte_zero_counters(prot_red);
                
                Monte_MC(prot_red);
                i_cycle++;
                
                /*
                this part is updated on 03/02/2007:
                sidechain->occ is now updated before written into the queue.
                */
                //cout << "after MC\n";
                for (ic=0;ic<prot_red.sidechain.size();ic++) {
                    float occ = (float) prot_red.sidechain[ic]->counter_accept / (float) (prot.param.monte_niter_cycle);
                    prot_red.sidechain[ic]->occ_table.push_back(occ);
                    prot_red.sidechain[ic]->occ = average(prot_red.sidechain[ic]->occ_table);
                    
                    //cout << prot_red.sidechain[ic]->uniqID << " " << prot_red.sidechain[ic]->occ << endl;
                    //cout << "ic=" << ic << "/" << prot_red.sidechain.size() <<endl;
                    /* back up occ at the end of last cycle*/
                    //prot_red.sidechain[ic]->occ_old = prot_red.sidechain[ic]->occ;
                    /* calc. avg occ of all finished cycles */
                    
                    /* back calculate avg occupancy within this cycle */
                    //float occ_of_this_cycle = prot_red.sidechain[ic]->occ * i_cycle - prot_red.sidechain[ic]->occ_old * (i_cycle-1);
                }
                
                if (i_red == 0) {
                    //Monte_group_conf_type(prot_red);
                    Monte_update_energy(prot_red); /* updating entropy term */
                    
                    /* update prot_red.E_state based on res.conf_w */
                    prot_red.E_state = prot_red.E_base;
                    for (int i_res=0;i_res<prot_red.res.size();i_res++) {
                        prot_red.E_state += prot_red.res[i_res].conf_w->E_self;
                        int ic = prot_red.res[i_res].conf_w->i_sidechain_prot_orig;
                        for (int j_res=i_res+1;j_res<prot_red.res.size();j_res++) {
                            int jc = prot_red.res[j_res].conf_w->i_sidechain_prot_orig;
                            prot_red.E_state += pairwise[ic][jc];
                        }
                    }
                }
                
                //cout << "after counting\n";
            }
            
            fp = fopen(MC_LOG,"a"); fprintf(fp,"   Actual  number of cycles to converge: %5d\n",i_cycle); fclose(fp);
            
            for (ic=0;ic<prot_w.sidechain.size();ic++) {
                prot_w.sidechain[ic]->counter_trial = 0; /* initialize because not all counters are updated later by prot_red */
            }
            
            /* upload occ and counter_trial of this reduced run to prot_w */
            jc = 0;
            for (ic=0;ic<prot_red.sidechain.size();ic++) {
                //prot_red.sidechain[ic]->occ = (float) prot_red.sidechain[ic]->counter_accept / (float) (i_cycle * prot.param.monte_niter_cycle);
                while (prot_red.sidechain[ic]->uniqID != prot_w.sidechain[jc]->uniqID) {
                    jc++;
                }
                prot_w.sidechain[jc]->occ = prot_red.sidechain[ic]->occ;
                //prot_w.sidechain[jc]->TS = prot_red.sidechain[ic]->TS;
                prot_w.sidechain[jc]->counter_trial = prot_red.sidechain[ic]->counter_trial;
            }
            
            for (ic=0;ic<prot_w.sidechain.size();ic++) {
                prot_w.sidechain[ic]->occ_table.push_back(prot_w.sidechain[ic]->occ);
                //prot_w.sidechain[ic]->TS_table.push_back(prot_w.sidechain[ic]->TS);
            }
            
            /* now update occ with avg value of all reduce runs */
            for (ic=0;ic<prot_w.sidechain.size();ic++) {
                prot_w.sidechain[ic]->occ = average(prot_w.sidechain[ic]->occ_table);
                //prot_w.sidechain[ic]->TS_table.push_back(prot_w.sidechain[ic]->TS);
            }
            
            //Monte_group_conf_type(prot_w);
            Monte_update_energy(prot_w);

            fp = fopen(MC_DETAIL,"a");
            fprintf(fp,"\n");
            fprintf(fp,"pH = %6.2f Eh = %6.2f, Unf. Energy = %8.2f Kcal/mol\n", ph, eh, prot.E_free_unfold);
            fprintf(fp,"pH = %6.2f Eh = %6.2f, Ave. Energy = %8.2f Kcal/mol\n", ph, eh, prot_red.E_accum / (float) (i_cycle * prot.param.monte_niter_cycle));
            fprintf(fp,"pH = %6.2f Eh = %6.2f, Min. Energy = %8.2f Kcal/mol\n", ph, eh, prot_red.E_min);
            fclose(fp);


            /*
            if (prot.param.titra_type == PH) titra_pnt = ph;
            else titra_pnt = eh;
            for (ic=0;ic<prot_w.sidechain.size();ic++) {
                fprintf(fp,"%6.2f %s occ=%5.3f trial=%8d/%8d %c\n",
                titra_pnt,
                prot_w.sidechain[ic]->uniqID,
                prot_w.sidechain[ic]->occ,
                prot_w.sidechain[ic]->counter_trial,
                i_cycle * prot.param.monte_niter_cycle,
                prot_w.sidechain[ic]->toggle);
            }
            */
            
            /*
            for (i_res=0;i_res<prot_red.res.size();i_res++) {
                if (prot_red.res[i_res].ngh.size()) prot_red.res[i_res].ngh.size();
            }
            if (prot_red.sidechain.size()) prot_red.sidechain.size();
            */
            i_red++;
        }
        
        /* update occ to prot */
        jc = 0;
        for (ic=0;ic<prot_w.sidechain.size();ic++) {
            while (prot_w.sidechain[ic]->uniqID != prot.sidechain[jc]->uniqID) jc++;
            prot.sidechain[jc]->occ = average(prot_w.sidechain[ic]->occ_table);
        }
        for (ic=0;ic<prot.sidechain.size();ic++) {
            prot.sidechain[ic]->occ_table.push_back(prot.sidechain[ic]->occ);
        }
        
        //Monte_group_conf_type(prot);
        Monte_update_energy(prot);
        
        for (ic=0;ic<prot.sidechain.size();ic++) {
            prot.sidechain[ic]->TS_table.push_back(prot.sidechain[ic]->TS);
        }

        timer_end = time(NULL);
        fp = fopen(MC_LOG,"a"); fprintf(fp,"Monte Carlo running time: %ld seconds.\n", timer_end-timer_start); fclose(fp);
        Monte_out(prot);
    }
    
    return 0;
}

int Monte_load_pairwise(Prot& prot)
{
    int   ic, jc;
    string  file_name;
    string  sbuff;
    string  stemp;
    int   i_res,i_conf,j_res,j_conf;
    
    /* declare memory */
    pairwise.resize(prot.sidechain.size());
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        pairwise[ic].resize(prot.sidechain.size(),0.0);
    }
    
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        if (prot.sidechain[ic]->uniqID.substr(3,2) == "DM") {
            continue;
        }
            
        string file_name = ENERGIES;
        file_name += "/" + prot.sidechain[ic]->uniqID + ".opp";
        load_opp(pairwise[ic], file_name.c_str(), prot, ic);
    }
    
    printf("   WARNING: Big difference (> %.1f%%) is reported\n", prot.param.monte_warn_pairwise);
    for (i_res=0; i_res<prot.res.size(); i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].conf.size(); i_conf++) {
            ic = prot.res[i_res].conf[i_conf].i_sidechain_prot;

            for (j_res=i_res; j_res<prot.res.size(); j_res++) {
                for (j_conf=1; j_conf<prot.res[j_res].conf.size(); j_conf++) {
                    jc = prot.res[j_res].conf[j_conf].i_sidechain_prot;
                    
                    if (i_res != j_res) {
                        /* warning for asymmetry */
                        float diff = fabs(pairwise[ic][jc]-pairwise[jc][ic]);
                        float avg  = (pairwise[ic][jc]+pairwise[jc][ic]);
                        if (diff*100./avg > prot.param.monte_warn_pairwise) {
                            if ( diff > 0.5 ) {
                                // cout << i_res << " " << j_res << endl;
                                printf("%s -> %s = %10.3f, <- = %10.3f\n",
                                prot.sidechain[ic]->uniqID.c_str(),
                                prot.sidechain[jc]->uniqID.c_str(),
                                pairwise[ic][jc], pairwise[jc][ic]);
                            }
                        }
                        
                        if (prot.param.monte_average_pairwise) {   /* average */
                            pairwise[ic][jc] = pairwise[jc][ic] = 0.5*(pairwise[ic][jc]+pairwise[jc][ic]);
                        }
                        else {      /* get the smaller */
                            pairwise[ic][jc] = pairwise[jc][ic] = fabs(pairwise[ic][jc]) > fabs(pairwise[jc][ic]) ? pairwise[jc][ic] : pairwise[ic][jc];
                        }
                    }
                    else {
                        pairwise[ic][jc] = pairwise[jc][ic] = 0.0; /* Necessary for fast E updating */
                    }
                }
            }
        }
    }
    
    return 0;
}

void Monte_get_biglist(Prot& prot)
{
    int i_res, j_res, i_conf, j_conf, ic, jc;
    bool added;
    
    for (i_res=0;i_res<prot.res.size();++i_res) {
        prot.res[i_res].ngh.clear();
        
        for (j_res=0;j_res<prot.res.size();j_res++) {
            if (i_res == j_res) continue;
            
            added = false;
            for (i_conf=1; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
                ic = prot.res[i_res].conf[i_conf].i_sidechain_prot_orig;
                for (j_conf=1;j_conf<prot.res[j_res].conf.size();j_conf++) {
                    jc = prot.res[j_res].conf[j_conf].i_sidechain_prot_orig;
                    
                    if (fabs(pairwise[ic][jc])>prot.param.big_pairwise) {
                        prot.res[i_res].ngh.push_back(&prot.res[j_res]);
                        added = true;
                        break;
                    }
                }
                if (added) break;
            }
        }
        
        prot.res[i_res].n_flip_max = (prot.res[i_res].ngh.size()+1) < prot.param.monte_flips ? (prot.res[i_res].ngh.size()+1) : prot.param.monte_flips;
    }
}

void Monte_check_toggle(Prot& prot)
{
    int i_res, i_conf, n_fixed;
    float fixed_occ;
    
    for (i_res=0;i_res<prot.res.size();i_res++) {
        fixed_occ = 0.;
        n_fixed = 0;
        for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
            if (prot.res[i_res].conf[i_conf].toggle == 'f') {
                n_fixed++;
                fixed_occ += prot.res[i_res].conf[i_conf].occ;
            }
        }
        if (fixed_occ > 1.) printf(" Warning! Total fixed occupancy on %s is over 1.\n",
            prot.res[i_res].resID().c_str());fflush(stdout);
        
        if (prot.res[i_res].conf.size() - n_fixed == 2) {
            for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
                if (prot.res[i_res].conf[i_conf].toggle != 'f') {
                    prot.res[i_res].conf[i_conf].toggle = 'f';
                    prot.res[i_res].conf[i_conf].occ = 1.-fixed_occ;
                    if (prot.res[i_res].conf[i_conf].occ < 0.) prot.res[i_res].conf[i_conf].occ = 0.;
                    n_fixed++;
                    fixed_occ += prot.res[i_res].conf[i_conf].occ;
                    break;
                }
            }
        }
        
        /* if total fixed occupancy is 1 fixed all the other free conformers */
        prot.res[i_res].fixed_occ = fixed_occ;
        if (prot.res[i_res].fixed_occ > (1.- 1e-6)) {
            if ((prot.res[i_res].conf.size() - n_fixed) > 1) {
                for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
                    if (prot.res[i_res].conf[i_conf].toggle != 'f') {
                        prot.res[i_res].conf[i_conf].toggle = 'f';
                        prot.res[i_res].conf[i_conf].occ = 0.;
                        n_fixed++;
                    }
                }
            }
        }
        
        /* if all confomers are fixed check if total fixed occupany is 1. */
        if (prot.res[i_res].conf.size() - n_fixed == 1) {
            if (fixed_occ < 1.) printf(" Warning! Total fixed occupancy on %s is smaller than 1.\n",
                prot.res[i_res].resID().c_str());fflush(stdout);
        }
    }
    return;
}

void Monte_set_toggle(Prot& prot)
{
    /* setup toggle based on MC occupancy */
    int i_res, i_conf;
    int n_fixed;
    
    for (i_res=0;i_res<prot.res.size();i_res++) {
        n_fixed = 0;
        for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
            if (prot.res[i_res].conf[i_conf].toggle == 'a') continue;
            if (prot.res[i_res].conf[i_conf].occ < prot.param.monte_reduce) {
                prot.res[i_res].conf[i_conf].toggle = 'f';
                prot.res[i_res].conf[i_conf].occ = 0.;
                n_fixed++;
            }
        }
        /* if only one conformer is free, fix it too */
        if ((prot.res[i_res].conf.size() - n_fixed) == 2) {
            for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
                if (prot.res[i_res].conf[i_conf].toggle != 'f') {
                    prot.res[i_res].conf[i_conf].toggle = 'f';
                    prot.res[i_res].conf[i_conf].occ = 1. - prot.res[i_res].fixed_occ;
                    prot.res[i_res].fixed_occ = 1.;
                    n_fixed++;
                    break;
                }
            }
        }
    }
}

Prot Monte_reduce(const Prot& prot)
{
    /* Uses toggle in prot to creat prot_red, reduced version of prot */
    Prot prot_red;
    int i_res, i_conf, j_res,j_conf;
    int ic, jc;
    
    prot_red = prot;
    
    for (i_res=prot_red.res.size()-1; i_res>=0; i_res--) {
        for (i_conf=prot_red.res[i_res].conf.size()-1; i_conf>=1; i_conf--) {
            if (prot_red.res[i_res].conf[i_conf].toggle == 'f') {       /* turned off */
                ic = prot_red.res[i_res].conf[i_conf].i_sidechain_prot_orig;
                
                prot_red.res[i_res].fixed_occ += prot_red.res[i_res].conf[i_conf].occ;
                
                /* Fix energy */
                prot_red.E_base += prot_red.res[i_res].conf[i_conf].E_self * prot_red.res[i_res].conf[i_conf].occ;
                
                for (j_res=0;j_res<prot_red.res.size();j_res++) {
                    if (i_res==j_res) continue;
                    for (j_conf=1;j_conf<prot_red.res[j_res].conf.size();j_conf++) {
                        jc = prot_red.res[j_res].conf[j_conf].i_sidechain_prot_orig;
                        prot_red.res[j_res].conf[j_conf].E_base += pairwise[ic][jc] * prot_red.res[i_res].conf[i_conf].occ;
                        prot_red.res[j_res].conf[j_conf].E_self += pairwise[ic][jc] * prot_red.res[i_res].conf[i_conf].occ;
                    }
                }
                
                /* delete this conformer */
                prot_red.res[i_res].del_conf(i_conf);
            }
        }
        
        /* delete the residue if no free conformer left */
        if (prot_red.res[i_res].conf.size() == 1) {
            prot_red.del_res(i_res);
        }
    }
    
    Monte_group_conf_type(prot_red);
    
    return prot_red;
}

void Monte_MC(Prot& prot)
{
    /* one cycle of monte carlo runs.
    All the counters are accumulated.
    03/02/2007: conf.occ not updated in this subroutine.
    Although conf.counter_accept is updated here, counter_iter is not kept track of */
    int k_weighted, k_res, k_conf, k_ngh;
    int n_flip, i_flip, i_iter;
    int ic_old, ic_new, jc_w;
    int j_res;
    double E_delta;
    
    i_iter = prot.param.monte_niter_cycle;
    while (i_iter) {
        
        /* select one residue */
        k_weighted = (int) (ran2(&dummy_num) * (float)prot.weighted_res_list.size());
        k_res = prot.weighted_res_list[k_weighted];
        //k_res = (int) (ran2(&dummy_num) * (double)prot.res.size());

        
        /* decide how many flips */
        n_flip = (int) ((ran2(&dummy_num) * (float)prot.res[k_res].n_flip_max)) + 1;
        
        /* get the list of residue to flip */
        flip_res[0] = &prot.res[k_res];
        for (i_flip=1;i_flip<n_flip;i_flip++) {
            int n_ngh = prot.res[k_res].ngh.size() - (i_flip-1);
            k_ngh = (int) (ran2(&dummy_num) * (float)n_ngh);
            flip_res[i_flip] = prot.res[k_res].ngh[k_ngh];
            swap(prot.res[k_res].ngh[k_ngh], prot.res[k_res].ngh[n_ngh-1]);
        }
        
        /* each residue pick new conformer */
        for (i_flip=0;i_flip<n_flip;i_flip++) {
            flip_res[i_flip]->conf_old = flip_res[i_flip]->conf_w;
            
            k_conf = (int)(ran2(&dummy_num) * (float)(flip_res[i_flip]->conf.size()-1)) + 1;
            flip_res[i_flip]->conf_new = &flip_res[i_flip]->conf[k_conf];
            flip_res[i_flip]->counter_trial++;
            flip_res[i_flip]->conf_new->counter_trial++;
        }
        
        /* calc. deltaE */
        E_delta = 0.;
        for (i_flip=0;i_flip<n_flip;i_flip++) {
            ic_old = flip_res[i_flip]->conf_old->i_sidechain_prot_orig;
            ic_new = flip_res[i_flip]->conf_new->i_sidechain_prot_orig;
            
            E_delta += (flip_res[i_flip]->conf_new->E_self - flip_res[i_flip]->conf_old->E_self);
            for (j_res=0; j_res<prot.res.size(); j_res++) {
                jc_w = prot.res[j_res].conf_w->i_sidechain_prot_orig;
                E_delta += (pairwise[ic_new][jc_w] - pairwise[ic_old][jc_w]);
            }
            
            /* real flip the active conformer */
            flip_res[i_flip]->conf_w = flip_res[i_flip]->conf_new;
        }
        
        /* Metropolis */
        if (exp(-beta*E_delta) > ran2(&dummy_num)) {
            prot.E_state += E_delta;
        }
        else {
            for (i_flip=0;i_flip<n_flip;i_flip++) {
                flip_res[i_flip]->conf_w = flip_res[i_flip]->conf_old;
            }
        }
        
        for (k_res=0;k_res<prot.res.size();k_res++) {
            prot.res[k_res].conf_w->counter_accept++;
        }
        if (prot.E_state < prot.E_min) prot.E_min = prot.E_state;
        prot.E_accum += prot.E_state;
        
        /*
        if (flip_res[0]->resSeq == 259) {
            if (flip_res[0]->conf_old->i_conf_res == 1) {
                if (flip_res[0]->conf_new->i_conf_res == 2) {
                    cout << "flip" << "\t";
                    cout << flip_res[0]->conf_old->uniqID << "\t";
                    cout << flip_res[0]->conf_old->E_self << "\t";
                    cout << flip_res[0]->conf_new->uniqID   << "\t";
                    cout << flip_res[0]->conf_new->E_self   << "\t";
                    cout << flip_res[0]->conf_w->uniqID   << "\t";
                    cout << flip_res[0]->conf_w->E_self   << "\t";
                    cout << E_delta << endl;
                }
            }
        }
        */
        //if (prot.param.monte_do_energy) do_free_energy(prot_p);
        i_iter--;
    }
    
    double E_chk = prot.E_base;
    for (int i_res=0;i_res<prot.res.size();i_res++) {
        E_chk += prot.res[i_res].conf_w->E_self;
        int ic = prot.res[i_res].conf_w->i_sidechain_prot_orig;
        //printf("checking, %s %d\n", prot.res[i_res].conf_w->uniqID.c_str(), ic);
        for (j_res=i_res+1;j_res<prot.res.size();j_res++) {
            int jc = prot.res[j_res].conf_w->i_sidechain_prot_orig;
            E_chk += pairwise[ic][jc];
        }
    }
    //printf("E_state=%10.4f,E_chk=%10.4f,diff=%.3e\n",prot.E_state,E_chk,prot.E_state-E_chk);

}

int Monte_out(Prot& prot)
{
    FILE *fp;
    int  ic, i_titra, i_res, i_conf;

    /* writing occ table */
    if (!(fp = fopen(OCC_TABLE, "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", OCC_TABLE);fflush(stdout);
        return USERERR;
    }
    if (prot.param.titra_type == PH) {
        fprintf(fp," ph           %s\n",vec2string(prot.titra_point,"%5.1f").c_str());
    }
    else {
        fprintf(fp," eh           %s\n",vec2string(prot.titra_point,"%5.0f").c_str());
    }
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        fprintf(fp,"%s %s\n",prot.sidechain[ic]->uniqID.c_str(), vec2string(prot.sidechain[ic]->occ_table, "%5.3f").c_str());
    }
    fclose(fp);

    /* writing entropy table */
    //Monte_update_TS(prot);
    if (!(fp = fopen("entropy.out", "w"))) {
        printf("   FATAL: Can not write occupancy to file \"entropy.out\"\n");fflush(stdout);
        return USERERR;
    }
    if (prot.param.titra_type == PH) {
        fprintf(fp," ph           %s\n",vec2string(prot.titra_point,"%5.1f").c_str());
    }
    else {
        fprintf(fp," eh           %s\n",vec2string(prot.titra_point,"%5.0f").c_str());
    }
    for (ic=0; ic<prot.sidechain.size(); ic++) {
        fprintf(fp,"%s %s\n",prot.sidechain[ic]->uniqID.c_str(), vec2string(prot.sidechain[ic]->TS_table, "%5.3f").c_str());
    }
    fclose(fp);
    
    /* writing sum_crg */
    prot.make_sum_crg();
    if (!(fp = fopen(TOT_CRG, "w"))) {
        printf("   FATAL: Can not write sumcrg to file \"%s\"\n", TOT_CRG);fflush(stdout);
        return USERERR;
    }
    if (prot.param.titra_type == PH) {
        fprintf(fp, " ph        ");
        fprintf(fp,"%s\n",vec2string(prot.titra_point,"%5.1f").c_str());
    }
    else {
        fprintf(fp, " eh        ");
        fprintf(fp,"%s\n",vec2string(prot.titra_point,"%5.0f").c_str());
    }
    
    for (i_res=0;i_res<prot.res.size();i_res++) {
        if (!prot.res[i_res].print_sumcrg) continue;
        fprintf(fp,"%s %s\n",prot.res[i_res].resID().c_str(),vec2string(prot.res[i_res].sum_crg,"%5.2f").c_str());
    }
    fprintf(fp,"Net_Charge %s\n",vec2string(prot.sum_crg,"%5.1f").c_str());
    fprintf(fp,"Protons    %s\n",vec2string(prot.sum_H,"%5.1f").c_str());
    fprintf(fp,"Electons   %s\n",vec2string(prot.sum_e,"%5.1f").c_str());
    fclose(fp);
    
    /* writing binding occupancy for water etc. */
    if (!(fp = fopen(TOT_OCC, "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", TOT_OCC);fflush(stdout);
        return USERERR;
    }
    if (prot.param.titra_type == PH) {
        fprintf(fp, " ph        ");
    }
    else {
        fprintf(fp, " eh        ");
    }
    fprintf(fp,"%s\n",vec2string(prot.titra_point,"%5.1f").c_str());
    
    for (i_res=0;i_res<prot.res.size();i_res++) {
        for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
            if (prot.res[i_res].conf[i_conf].uniqID.substr(3,2) == "DM") break;
        }
        if (i_conf >= prot.res[i_res].conf.size()) continue;

        fprintf(fp, "%s", prot.res[i_res].resID().c_str());
        for (i_titra=0; i_titra<prot.res[i_res].conf[1].occ_table.size(); i_titra++) {
            float occ = 1.;
            for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
                if (prot.res[i_res].conf[i_conf].uniqID.substr(3,2) == "DM")
                    occ -= prot.res[i_res].conf[i_conf].occ_table[i_titra];
            }
            fprintf(fp, " %5.2f", occ);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    
    if (!(fp = fopen(CURVE_FITTING, "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", CURVE_FITTING);fflush(stdout);
        return USERERR;
    }
    if (prot.param.titra_type == PH) {
        fprintf(fp, "Residue            pK      n(slope)  1000*chi2\n");
    }
    else {
        fprintf(fp, "Residue            Em      n(slope)  1000*chi2\n");
    }
    for (i_res=0;i_res<prot.res.size();i_res++) {
        if (!prot.res[i_res].print_sumcrg) continue;
        Pka pka;
        pka = curve_fitting(prot.titra_point, prot.res[i_res].sum_crg);
        if (pka.valid) {
            if (prot.param.titra_type == EH) pka.n *= 58.;
            fprintf(fp,"%s    %9.3f %9.3f %9.3f",prot.res[i_res].resID().c_str(), pka.pK, pka.n, 1000.*pka.chi2);
            if (pka.out_of_range || fabs(pka.n) > 1.1) fprintf(fp,"   Warning! Titration is not complete, check sumcrg.");
            else if (fabs(pka.n) < 0.8 || pka.chi2 > 2e-3) fprintf(fp,"   Warning! This titration may be coupled to other groups, check sumcrg.");
            fprintf(fp,"\n");
        }
        else {
            fprintf(fp,"%s         Titration is out of range.\n",prot.res[i_res].resID().c_str());
        }
    }
    fclose(fp);
    return 0;
}

void Monte_zero_counters(Prot& prot)
{
    int i_res,i_conf;
    prot.E_accum = 0.;
    for (i_res=0;i_res<prot.res.size();i_res++) {
        prot.res[i_res].counter_trial = 0;
        for (i_conf=1;i_conf<prot.res[i_res].conf.size();i_conf++) {
            prot.res[i_res].conf[i_conf].counter_trial = 0;
            prot.res[i_res].conf[i_conf].counter_accept = 0;
        }
    }
}

void Monte_update_TS(Prot& prot)
{
    for (int i_res=0; i_res<prot.res.size(); ++i_res) {
        for (int i_conf = 1; i_conf<prot.res[i_res].conf.size(); ++i_conf) {
            prot.res[i_res].conf[i_conf].TS = 0.;
        }
        
        for (int i_conf_type=0; i_conf_type<prot.res[i_res].conf_type.size(); ++i_conf_type) {
            double sum_occ = 0.;  //total occupancy of this conf type, used for normalization 
            
            for (int i_slot=0;i_slot<prot.res[i_res].conf_type[i_conf_type].size();++i_slot) {
                int i_conf = prot.res[i_res].conf_type[i_conf_type][i_slot];
                sum_occ += prot.res[i_res].conf[i_conf].occ;
            }
            
            if (sum_occ > 1e-5) {
                double TS = 0.;
                for (int i_slot=0;i_slot<prot.res[i_res].conf_type[i_conf_type].size();++i_slot) {
                    int i_conf = prot.res[i_res].conf_type[i_conf_type][i_slot];
                    double renormalized_occ = prot.res[i_res].conf[i_conf].occ/sum_occ;
                    if (renormalized_occ > 1e-5) {
                        TS -= renormalized_occ*log(renormalized_occ) * KT2KCAL;
                    }
                }
                
                for (int i_slot=0;i_slot<prot.res[i_res].conf_type[i_conf_type].size();++i_slot) {
                    int i_conf = prot.res[i_res].conf_type[i_conf_type][i_slot];
                    prot.res[i_res].conf[i_conf].TS = TS;
                }
            }
        }
        
    }
}

void Monte_group_conf_type(Prot& prot)
{
    for (int i_res=0; i_res<prot.res.size(); ++i_res) {
        prot.res[i_res].conf_type.clear();
        for (int i_conf=1;i_conf<prot.res[i_res].conf.size();++i_conf) {
            int i_conf_type;
            for (i_conf_type=0; i_conf_type<prot.res[i_res].conf_type.size(); ++i_conf_type) {
                int k_conf = prot.res[i_res].conf_type[i_conf_type][0];
                if (prot.res[i_res].conf[i_conf].e != prot.res[i_res].conf[k_conf].e) continue;
                if (prot.res[i_res].conf[i_conf].H != prot.res[i_res].conf[k_conf].H) continue;
                if ((prot.res[i_res].conf[i_conf].confName.substr(3,2) == "DM") !=
                    (prot.res[i_res].conf[k_conf].confName.substr(3,2) == "DM")) continue; /* either both are dummy, or neither */
                /* add to the same conformer type */
                prot.res[i_res].conf_type[i_conf_type].push_back(i_conf);
                break;
            }
            if (i_conf_type >= prot.res[i_res].conf_type.size()) {
                vector<int> buff_vec(1);
                buff_vec[0] = i_conf;
                /* create a new conformer type */
                prot.res[i_res].conf_type.push_back(buff_vec);
            }
        }
    }
}

void Monte_update_energy(Prot& prot)
{   /* update entropy based on conf.occ, then update conf.E_self */
    Monte_group_conf_type(prot);
    Monte_update_TS(prot);
    
    /* calc. pH and Eh dependent self-energy of each conformer */
    float ph = prot.param.titra_ph0;
    float eh = prot.param.titra_eh0;
    if (prot.param.titra_type == PH) {
        ph = prot.titra_point[prot.titra_point.size()-1];
    }
    else {
        eh = prot.titra_point[prot.titra_point.size()-1];
    }        
    for (int ic=0; ic<prot.sidechain.size(); ic++) {
        prot.sidechain[ic]->E_ph =
        prot.param.monte_temp/ROOMT * prot.sidechain[ic]->H * (ph-prot.sidechain[ic]->pKa) * PH2KCAL;
        prot.sidechain[ic]->E_eh =
        prot.param.monte_temp/ROOMT * prot.sidechain[ic]->e * (eh-prot.sidechain[ic]->Em) * MV2KCAL;
        
        prot.sidechain[ic]->E_self
        = prot.sidechain[ic]->E_base
        + prot.sidechain[ic]->E_vdw0
        + prot.sidechain[ic]->E_vdw1
        + prot.sidechain[ic]->E_epol
        + prot.sidechain[ic]->E_tors
        + prot.sidechain[ic]->E_dsolv
        + prot.sidechain[ic]->E_extra
        + prot.sidechain[ic]->E_ph
        + prot.sidechain[ic]->E_eh
        + prot.sidechain[ic]->TS;       /* entropy correction */
        //printf("update ic=%3d, %s, %8.3f\n",ic, prot.sidechain[ic]->uniqID.c_str(), prot.sidechain[ic]->E_self);
    }
    
    /* update prot.E_state based on res.conf_w 
    prot.E_state = prot.E_base;
    for (int i_res=0;i_res<prot.res.size();i_res++) {
        prot.E_state += prot.res[i_res].conf_w->E_self;
        int ic = prot.res[i_res].conf_w->i_sidechain_prot_orig;
        for (int j_res=i_res+1;j_res<prot.res.size();j_res++) {
            int jc = prot.res[j_res].conf_w->i_sidechain_prot_orig;
            prot.E_state += pairwise[ic][jc];
        }
    }
    */
}

