#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mcce.hpp"

float get_chi2(float param[]);
float get_chi2_2pK(float param[]);
float get_chi2_N_pK(float param[]);

float  _tmp_sumcrg_start;
float  _tmp_sumcrg_end;
int    _tmp_n_titra;
float* _tmp_titra_point;
float* _tmp_sumcrg;
int    _N_pK;

Pka curve_fitting(vector <float> titra_point, vector <float> sumcrg_table)
{
    return curve_fitting(titra_point, sumcrg_table, 0.0, 0.0);
}

Pka curve_fitting(vector <float> titra_point, vector <float> sumcrg_table, float shift_pK, float shift_n)
{
    /* initialize the simplex and set up optimization */
    float** param;
    float chi2[3];
    int   i, n_iter, i_titra;
    Pka   pka;
    float same_occ_thr = 0.05;
    
    pka.valid = false;
    if (titra_point.size() != sumcrg_table.size()) {
        return pka;
    }
    
    float max_occ = max(sumcrg_table);
    float min_occ = min(sumcrg_table);
    
    if (max_occ - min_occ < 0.1) return pka;
    
    _tmp_n_titra     = titra_point.size();
    _tmp_titra_point = (float *)malloc(_tmp_n_titra*sizeof(float));
    _tmp_sumcrg      = (float *)malloc(_tmp_n_titra*sizeof(float));
    for (i_titra=0;i_titra<_tmp_n_titra; ++i_titra) {
        _tmp_titra_point[i_titra] = titra_point[i_titra];
        _tmp_sumcrg[i_titra]      = sumcrg_table[i_titra];
    }
    
    /* try to guess pK and n */
    float hlf_occ = (max_occ + min_occ) /2.;
    i_titra = nearest_index(hlf_occ, sumcrg_table);
    float guessed_pK = titra_point[i_titra];
    if (average(subvec(sumcrg_table, 0, i_titra))
    < average(subvec(sumcrg_table, i_titra, sumcrg_table.size()))) {
        _tmp_sumcrg_start = min_occ;
        _tmp_sumcrg_end = max_occ;
    }
    else {
        _tmp_sumcrg_start = max_occ;
        _tmp_sumcrg_end = min_occ;
    }
    float titra_begin = titra_point[rfind_match_index(_tmp_sumcrg_start, sumcrg_table, same_occ_thr)];
    float titra_end = titra_point[find_match_index(_tmp_sumcrg_end, sumcrg_table, same_occ_thr)];
    if (fabs(titra_end - titra_begin) < 1e-8) return pka;
    
    float guessed_n  = 4.0/(titra_end - titra_begin);
    if (fabs(shift_pK) < 1e-8) shift_pK = (titra_end - titra_begin)/4.0;
    if (fabs(shift_n) < 1e-8) shift_n = -guessed_n/2.0;
    
    /* check if titration is complete */
    pka.out_of_range = false;
    if ( fabs(_tmp_sumcrg_start - rintf(_tmp_sumcrg_start) ) > same_occ_thr ) {
        if (fabs(_tmp_sumcrg_start-sumcrg_table[0]) < same_occ_thr) {
            pka.out_of_range = true;
        }
    }
    if ( fabs(_tmp_sumcrg_end - rintf(_tmp_sumcrg_end) ) > same_occ_thr ) {
        if (fabs(_tmp_sumcrg_end-sumcrg_table[_tmp_n_titra-1]) < same_occ_thr) {
            pka.out_of_range = true;
        }
    }
    
    param = (float**)malloc(3 * sizeof(float *));
    for (i=0;i<3;i++) param[i] = (float *)malloc(2 * sizeof(float));
    
    param[0][0] = guessed_pK;           param[0][1] = guessed_n;        chi2[0] = get_chi2(param[0]);
    param[1][0] = guessed_pK+shift_pK;  param[1][1] = guessed_n;        chi2[1] = get_chi2(param[1]);
    param[2][0] = guessed_pK;           param[2][1] = guessed_n+shift_n;chi2[2] = get_chi2(param[2]);
    
    dhill(param, chi2, 2, 1e-9, &get_chi2, &n_iter);
    
    /*  fit this eq. sumcrg = sumstart + (sumend - sumstart) * exp(n(pH-pK))/(1+exp(n(pH-pK)))  */
    
    pka.pK =   param[0][0];
    pka.n  =   param[0][1];
    pka.chi2 = chi2[0];
    
    for (i=0;i<3;i++)
        free(param[i]);
    free(param);
    
    if ( (min(sumcrg_table)+max(sumcrg_table)) < 0.) {
        pka.n = -fabs(pka.n);
    }
    else {
        pka.n = fabs(pka.n);
    }
    
    pka.valid = true;
    return pka;
}

float get_chi2(float param[])
{
    /* param[0] = pKa, param[1] = n */
    float   chi2 = 0.0;
    double  y_fit;
    double  T;
    int     i;
    
    for (i=0; i<_tmp_n_titra; i++) {
        T = 2.303*param[1]*(_tmp_titra_point[i]-param[0]);
        
        if (T<50.)
            y_fit = _tmp_sumcrg_start + (_tmp_sumcrg_end - _tmp_sumcrg_start) * exp(T)/(1.0+exp(T));
        else
            y_fit = _tmp_sumcrg_end;
        
        chi2 += (y_fit-_tmp_sumcrg[i])*(y_fit-_tmp_sumcrg[i]);
    }
    
    return chi2;
}

vector<Pka> curve_fitting_2pK(vector <float> titra_point, vector <float> sumcrg_table)
{
    /* initialize the simplex and set up optimization */
    float** param;
    float chi2[6];
    int   i, n_iter, i_titra;
    vector<Pka>   pka(2);
    //float same_occ_thr = 0.05;
    
    pka[0].valid = false;
    pka[1].valid = false;
    if (titra_point.size() != sumcrg_table.size()) {
        return pka;
    }
    
    float max_occ = max(sumcrg_table);
    float min_occ = min(sumcrg_table);
    
    if (max_occ - min_occ < 0.1) return pka;
    
    _tmp_n_titra     = titra_point.size();
    _tmp_titra_point = (float *)malloc(_tmp_n_titra*sizeof(float));
    _tmp_sumcrg      = (float *)malloc(_tmp_n_titra*sizeof(float));
    for (i_titra=0;i_titra<_tmp_n_titra; ++i_titra) {
        _tmp_titra_point[i_titra] = titra_point[i_titra];
        _tmp_sumcrg[i_titra]      = sumcrg_table[i_titra];
    }
    
    /* try to guess pK and n */
    float hlf_occ = (max_occ + min_occ) /2.;
    i_titra = nearest_index(hlf_occ, sumcrg_table);
    //float guessed_pK1 = titra_point[i_titra]-1.;
    //float guessed_pK2 = titra_point[i_titra]+1.;
    
    /* decide if titration is going up or down */
    if (average(subvec(sumcrg_table, 0, i_titra))
    < average(subvec(sumcrg_table, i_titra, sumcrg_table.size()))) {
        _tmp_sumcrg_start = min_occ;
        _tmp_sumcrg_end = max_occ;
    }
    else {
        _tmp_sumcrg_start = max_occ;
        _tmp_sumcrg_end = min_occ;
    }
    
    /* check if titration is too rapid */
    //float titra_begin = titra_point[rfind_match_index(_tmp_sumcrg_start, sumcrg_table, same_occ_thr)];
    //float titra_end = titra_point[find_match_index(_tmp_sumcrg_end, sumcrg_table, same_occ_thr)];
    //if (fabs(titra_end - titra_begin) < 1e-8) return pka;
    
    //float guessed_n  = 1.0;
    //float shift_pK = 0.1;
    //float shift_n = 0.1;
    
    /* check if titration is complete 
    pka[0].out_of_range = false;
    if ( fabs(_tmp_sumcrg_start - rintf(_tmp_sumcrg_start) ) > same_occ_thr ) {
        if (fabs(_tmp_sumcrg_start-sumcrg_table[0]) < same_occ_thr) {
            pka.out_of_range = true;
        }
    }
    if ( fabs(_tmp_sumcrg_end - rintf(_tmp_sumcrg_end) ) > same_occ_thr ) {
        if (fabs(_tmp_sumcrg_end-sumcrg_table[_tmp_n_titra-1]) < same_occ_thr) {
            pka.out_of_range = true;
        }
    }
    */
    
    param = (float**)malloc(6 * sizeof(float *));
    for (i=0;i<6;i++) param[i] = (float *)malloc(5 * sizeof(float));
    
    /* find closest set */
    chi2[1] = 99999.;
    //param[0][4] = 0.5;
    for (param[0][0] =  0.; param[0][0] < 14.1; param[0][0] +=0.5) {
    for (param[0][2] =  0.; param[0][2] < 14.1; param[0][2] +=0.5) {
    for (param[0][1] = -1.; param[0][1] <  1.1; param[0][1] +=0.2) {
    for (param[0][3] = -1.; param[0][3] <  1.1; param[0][3] +=0.2) {
    for (param[0][4] =  0.; param[0][4] <  1.1; param[0][4] +=0.2) {
        
        chi2[0] = get_chi2_2pK(param[0]);
        if (chi2[0] < chi2[1]) {
            chi2[1] = chi2[0];
            for (int i=0;i<5;i++) {
                param[1][i] = param[0][i];
            }
        }
    }
    }
    }
    }
    }
    
    for (int i=0;i<5;i++) {
    for (int j=0;j<6;j++) {
        param[j][i] = param[1][i];
    }
    }
    for (int i=0;i<5;i++) {
        param[i+1][i] += 0.01;
    }
    
    for (int i=0;i<6;i++) {
        chi2[i] = get_chi2_2pK(param[i]);
    }
    
    dhill(param, chi2, 5, 1e-9, &get_chi2_2pK, &n_iter);
    
    /*  fit this eq. sumcrg = sumstart + (sumend - sumstart) * exp(n(pH-pK))/(1+exp(n(pH-pK)))  */
    
    pka[0].pK =   param[0][0];
    pka[0].n  =   param[0][1];
    pka[0].coef = param[0][4];
    pka[0].chi2 = chi2[0];
    
    pka[1].pK =   param[0][2];
    pka[1].n  =   param[0][3];
    pka[1].coef = _tmp_sumcrg_end - _tmp_sumcrg_start - param[0][4];
    pka[1].chi2 = chi2[0];
    
    for (i=0;i<6;i++)
        free(param[i]);
    free(param);
    
    if ( (min(sumcrg_table)+max(sumcrg_table)) < 0.) {
        pka[0].n = -fabs(pka[0].n);
        pka[1].n = -fabs(pka[1].n);
    }
    else {
        pka[0].n = fabs(pka[0].n);
        pka[1].n = fabs(pka[1].n);
    }
    
    pka[0].valid = true;
    pka[1].valid = true;
    return pka;
}

float get_chi2_2pK(float param[])
{
    /* param[0] = pK1, param[1] = n1 */
    /* param[2] = pK2, param[3] = n2 */
    /* param[4] = coefficient1 */
    
    float   chi2 = 0.0;
    double  y_fit;
    double  T1,T2;
    int     i;
    
    for (i=0; i<_tmp_n_titra; i++) {
        T1 = 2.303*param[1]*(_tmp_titra_point[i]-param[0]);
        T2 = 2.303*param[3]*(_tmp_titra_point[i]-param[2]);
        
        y_fit = _tmp_sumcrg_start;
        if (T1<50.)
            y_fit += param[4] * exp(T1)/(1.0+exp(T1));
        else
            y_fit += param[4];
        
        if (T2<50.)
            y_fit += (_tmp_sumcrg_end - _tmp_sumcrg_start - param[4]) * exp(T2)/(1.0+exp(T2));
        else
            y_fit += (_tmp_sumcrg_end - _tmp_sumcrg_start - param[4]);
        
        chi2 += (y_fit-_tmp_sumcrg[i])*(y_fit-_tmp_sumcrg[i]);
    }
    
    return chi2;
}

vector<Pka> curve_fitting_N_pK(vector <float> titra_point, vector <float> sumcrg_table, int& N_pK)
{
    /* initialize the simplex and set up optimization */
    /* For N pKs, each has three parameter: */
    /* param[0] = coefficient; param[1] = pK; param[2] = n */
    
    float** param;
    float* chi2;
    int   n_iter, i_titra;
    vector<Pka>   pka;
    float same_occ_thr = 0.05;
    int N_param;
    N_param = 3*N_pK;
    _N_pK = N_pK;
    
    pka.resize(N_pK);
    for (int i=0;i<N_pK;i++) {
        pka[i].valid = false;
    }
    if (titra_point.size() != sumcrg_table.size()) {
        return pka;
    }
    
    float max_occ = max(sumcrg_table);
    float min_occ = min(sumcrg_table);
    
    if (max_occ - min_occ < 0.1) return pka;
    
    _tmp_n_titra     = titra_point.size();
    _tmp_titra_point = (float *)malloc(_tmp_n_titra*sizeof(float));
    _tmp_sumcrg      = (float *)malloc(_tmp_n_titra*sizeof(float));
    for (i_titra=0;i_titra<_tmp_n_titra; ++i_titra) {
        _tmp_titra_point[i_titra] = titra_point[i_titra];
        _tmp_sumcrg[i_titra]      = sumcrg_table[i_titra];
    }
    
    /* try to guess pK and n */
    float hlf_occ = (max_occ + min_occ) /2.;
    i_titra = nearest_index(hlf_occ, sumcrg_table);
    //float guessed_pK1 = titra_point[i_titra]-1.;
    //float guessed_pK2 = titra_point[i_titra]+1.;
    
    /* decide if titration is going up or down */
    if (average(subvec(sumcrg_table, 0, i_titra))
    < average(subvec(sumcrg_table, i_titra, sumcrg_table.size()))) {
        _tmp_sumcrg_start = min_occ;
        _tmp_sumcrg_end = max_occ;
    }
    else {
        _tmp_sumcrg_start = max_occ;
        _tmp_sumcrg_end = min_occ;
    }
    
    /* check if titration is too rapid for fitting */
    float titra_begin = titra_point[rfind_match_index(_tmp_sumcrg_start, sumcrg_table, same_occ_thr)];
    float titra_end = titra_point[find_match_index(_tmp_sumcrg_end, sumcrg_table, same_occ_thr)];
    if (fabs(titra_end - titra_begin) < 1e-8) return pka;
    
    /* check if titration is complete 
    pka[0].out_of_range = false;
    if ( fabs(_tmp_sumcrg_start - rintf(_tmp_sumcrg_start) ) > same_occ_thr ) {
        if (fabs(_tmp_sumcrg_start-sumcrg_table[0]) < same_occ_thr) {
            pka.out_of_range = true;
        }
    }
    if ( fabs(_tmp_sumcrg_end - rintf(_tmp_sumcrg_end) ) > same_occ_thr ) {
        if (fabs(_tmp_sumcrg_end-sumcrg_table[_tmp_n_titra-1]) < same_occ_thr) {
            pka.out_of_range = true;
        }
    }
    */
    
    chi2  = (float *)  malloc((N_param+1) * sizeof(float *));
    param = (float **) malloc((N_param+1) * sizeof(float *));
    for (int i=0;i<N_param+1;i++) param[i] = (float *)malloc(N_param * sizeof(float));

    /* find closest set */
    float min_chi2 = 99999.;
    float *bes_fit_param;
    bes_fit_param = (float *)malloc(N_param * sizeof(float));
    
    for (int i_pK=0; i_pK<N_pK; i_pK++) {
        for (param[0][i_pK*3] = -1.; param[0][i_pK*3] <  1.1; param[0][i_pK*3] +=0.2) {  /* coefficient of the i-th pK */
            for (param[0][i_pK*3+1] =  0.; param[0][i_pK*3+1] < 14.1; param[0][i_pK*3+1] +=0.5) {  /* pK */
                for (param[0][i_pK*3+2] = -1.; param[0][i_pK*3+2] <  1.1; param[0][i_pK*3+2] +=0.2) { /* n */
                    
                    chi2[0] = get_chi2_N_pK(param[0]);
                    if (chi2[0] < min_chi2) {
                        min_chi2 = chi2[0];
                        for (int i=0;i<N_param;i++) {
                            bes_fit_param[i] = param[0][i];
                        }
                    }
                }
            }
        }
    }
    
    for (int i=0;i<N_param+1;i++) {
    for (int j=0;j<N_param;j++) {
        param[i][j] = bes_fit_param[j];
    }
    }
    for (int i_pK=0;i_pK<N_pK;i_pK++) {
        param[i_pK*3+1][i_pK*3]   +=0.1;
        param[i_pK*3+2][i_pK*3+1] +=0.5;
        param[i_pK*3+3][i_pK*3+2] -=0.1;
    }
    
    for (int i=0;i<N_param+1;i++) {
        chi2[i] = get_chi2_N_pK(param[i]);
    }
    
    dhill(param, chi2, N_param, 1e-9, &get_chi2_N_pK, &n_iter);
    
    /*  fit this eq. sumcrg = sumstart + (sumend - sumstart) * exp(n(pH-pK))/(1+exp(n(pH-pK)))  */
    
    for (int i_pK=0;i_pK<N_pK;i_pK++) {
        pka[i_pK].coef = param[0][i_pK*3];
        pka[i_pK].pK =   param[0][i_pK*3+1];
        pka[i_pK].n  =   param[0][i_pK*3+2];
        pka[i_pK].chi2 = chi2[0];
    }
    
    for (int i=0;i<N_param+1;i++)
        free(param[i]);
    free(param);
    
    if ( (min(sumcrg_table)+max(sumcrg_table)) < 0.) {
        for (int i_pK=0;i_pK<N_pK;i_pK++) {
            pka[i_pK].n = -fabs(pka[i_pK].n);
        }
    }
    else {
        for (int i_pK=0;i_pK<N_pK;i_pK++) {
            pka[i_pK].n = fabs(pka[i_pK].n);
        }
    }
    
    pka[0].valid = true;
    pka[1].valid = true;
    return pka;
}

float get_chi2_N_pK(float param[])
{
    /* param[0] = pK1, param[1] = n1 */
    /* param[2] = pK2, param[3] = n2 */
    /* param[4] = coefficient1 */
    
    float   chi2 = 0.0;
    double  y_fit; /* y value of the fitting curve */
    int     i;
    
    for (i=0; i<_tmp_n_titra; i++) { /* loop over all titration points */
        
        y_fit = _tmp_sumcrg_start;
        for (int i_pK=0; i_pK<_N_pK; i_pK++) { /* loop over all fitted pK */
            float coefficient = param[3*i_pK];
            float pK  = param[3*i_pK+1];
            float n   = param[3*i_pK+2];
            float pH  = _tmp_titra_point[i];
            
            double T = 2.303*n*(pH-pK);

            if (T<50.)
                y_fit += coefficient * exp(T)/(1.0+exp(T));
            else
                y_fit += coefficient;
        }
        
        /* accumulate chi2 */
        chi2 += (y_fit-_tmp_sumcrg[i])*(y_fit-_tmp_sumcrg[i]);
    }
    
    return chi2;
}

#define TINY 1.0E-10
#define NMAX 50000
#define GET_PSUM for (j=0; j<ndim; j++) {\
                        for (sum=0.0, i=0; i<mpts; i++) sum += p[i][j];\
		                 psum[j] = sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk)
{
    float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac);
    int i, ihi, ilo, inhi,j, mpts = ndim+1;
    float rtol,sum,swap,ysave,ytry,*psum;
    
    psum = (float *)malloc(ndim * sizeof(float));
    *nfunk = 0;
    GET_PSUM
    
    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1,0):(inhi = 0,1);
        for(i=0; i<mpts;i++) {
            if(y[i] <= y[ilo]) ilo=i;
            if(y[i] > y[ihi]) {
                inhi = ihi;
                ihi  = i;
            }
            else if (y[i] > y[inhi] && i != ihi) inhi = i;
        }
        
        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        
        if(rtol < ftol || *nfunk >= NMAX) {
            SWAP(y[0],y[ilo])
            for (i=0; i<ndim; i++) SWAP(p[0][i],p[ilo][i])
                break;
        }
        
        *nfunk += 2;
        
        ytry = dhtry(p,y,psum,ndim,funk,ihi,-1.0);
        
        if (ytry <= y[ilo]) ytry=dhtry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry=dhtry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0; i<mpts; i++) {
                    if (i !=ilo) {
                        for (j=0; j<ndim; j++)
                            p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);
                        y[i] = (*funk)(psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        }
        else --(*nfunk);
    }
    free(psum);
}

float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac)
{
    int j;
    float fac1, fac2, ytry, *ptry;
    
    ptry = (float *)malloc(ndim*sizeof(float));
    fac1 = (1.0-fac)/ndim;
    fac2 = fac1 - fac;
    for (j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry = (*funk)(ptry);   /* evaluate the function at the trial point */
    if (ytry < y[ihi]) {    /* if it's better than the highest, then replace the highest */
        y[ihi] = ytry;
        for (j=0; j<ndim; j++) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    free(ptry);
    return ytry;
}

