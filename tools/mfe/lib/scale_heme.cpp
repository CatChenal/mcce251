int main(int argc, char** argv) {
    Prot prot;
    
    /* Load head3.lst */
    try {
        prot.load_conflist(FN_CONFLIST3);
    }
    catch (string msg) {    /* if there is error in loading conflist */
        return USERERR;
    }
    prot.load_pdb("step3_out.pdb");

    /* Load pairwise */
    scale_ele = 1.0;
    scale_vdw = 0.0;
    if (Monte_load_pairwise(prot)) return USERERR;
    
    /* rescale elec. pair for heme */




    if (prot.sidechain[kconf]->confName == "HAN+W"
        ||  prot.sidechain[kconf]->confName == "HA3+W"
    ||  prot.sidechain[ic]->confName == "HAN+W"
    ||  prot.sidechain[ic]->confName == "HA3+W" ) {
        
        /*
        if (prot.sidechain[kconf]->res->resName == "ARG" ||
            prot.sidechain[kconf]->res->resName == "HIS"
        ||  prot.sidechain[kconf]->res->resName == "CUB"
        ||  prot.sidechain[ic]->res->resName == "ARG"
        ||  prot.sidechain[ic]->res->resName == "HIS"
        ||  prot.sidechain[ic]->res->resName == "CUB" ) {
        }
        */
        
        if (ele_pair > 1) {
            
            cout << prot.sidechain[kconf]->uniqID << " <-> "
            << prot.sidechain[ic]->uniqID << " "
            << ele_pair << " --> ";
            
            float heme_scale = 1-0.4*ele_pair/10.;
            if (heme_scale < 0.4) heme_scale = 0.4;
            ele_pair = ele_pair * heme_scale;
            
            cout << ele_pair << endl;
        }
    }

