#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "mcce.hpp"

int main(int argc, char** argv) {
    Prot prot1, prot2;
    
    /* parsing command */
    if (argc < 3) {
        cout << "Usage:\n";
        cout << "   " << argv[0] << " file1 file2\n";
        cout << "   file1 is the head3.lst used to get flags and occupancies\n";
        cout << "   file2 is the head3.lst used to get conformer list and numbers\n";
        return -1;
    }
    
    prot1.load_conflist(argv[1]);
    prot1.index_on();
    
    prot2.load_conflist(argv[2]);
    prot2.index_on();
    
    for (int ic = 0; ic<prot1.sidechain.size(); ++ic) {
        for (int jc = 0; jc<prot2.sidechain.size(); ++jc) {
            if (prot1.sidechain[ic]->uniqID == prot2.sidechain[jc]->uniqID) {
                prot2.sidechain[jc]->toggle = prot1.sidechain[ic]->toggle;
                prot2.sidechain[jc]->occ   = prot1.sidechain[ic]->occ;
                break;
            }
        }
    }
    
    prot2.write_conflist("head3.lst.new");
    return 0;
}
