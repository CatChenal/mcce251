
int load_pairwise()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, ".")<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    /* scan conformers to see if all pw were calculated */
    for (i=0; i<ematrix.n; i++) {
        if (!ematrix.conf[i].on && (ematrix.conf[i].uniqID[3] != 'D' || ematrix.conf[i].uniqID[4] != 'M')) {
           printf("      Incompleted delphi run, the first place detected at %s\n ", ematrix.conf[i].uniqID);
           return USERERR;
        }
    }
    
    
    if (!(pairwise = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    for (kc=0; kc<ematrix.n; kc++) {
        if (!(pairwise[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           /* proprocessing */
           if ((ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw) > 999.0) {
              pairwise[i][j] = pairwise[j][i] = 999.0;
           }
           else {
              pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele \
                                                +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           }
        }
    }

    /* DEBUG print pairwise table 
    int kr;
    printf("                ");
    for (kc=0; kc<conflist.n_conf; kc++) {
        printf("%16s", conflist.conf[kc].uniqID);
    }
    printf("\n");
    for (kr=0; kr<conflist.n_conf; kr++) {
        printf("%s ", conflist.conf[kr].uniqID);
        for (kc=0; kc<conflist.n_conf; kc++) {
            printf("%16.3f", pairwise[kr][kc]);
        }
        printf("\n");
    }
    */

    
    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

