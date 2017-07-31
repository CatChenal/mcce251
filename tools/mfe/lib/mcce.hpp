/********************************
* MCCE library header file
********************************/
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "math3d.hpp"

using namespace std;
using namespace Math3D;

#define USERERR     -1

/* at 298.15K = 25C */
#define ROOMT        298.15
#define KCAL2KT      1.688
#define KT2KCAL      0.593
#define PH2KCAL      1.364
#define PH2MV        59.16
#define MV2KCAL      0.023

/* IO */
#define MC_LOG       "mc_out"
#define MC_DETAIL    "fort.36"
#define OCC_TABLE    "fort.38"
#define TOT_CRG      "sum_crg.out"
#define TOT_OCC      "sum_hoh.out"
#define CURVE_FITTING "pK.out"
#define FN_CONFLIST1 "head1.lst"
#define FN_CONFLIST2 "head2.lst"
#define FN_CONFLIST3 "head3.lst"
#define STEP1_OUT    "step1_out.pdb"
#define STEP2_OUT    "step2_out.pdb"
#define STEP3_OUT    "step3_out.pdb"
#define ENERGIES     "energies"
#define ENERGY_TABLE "energies.opp"

class Param;
class Prot;
class Res;
class Conf;
class Atom;

enum    HetType {ATOM,HETATM};
enum    TitraType {PH,EH};  // Need CH also?

bool sameID(Conf &conf1, Conf &conf2);

class Param_value {
    public:

    private:

	    string  _key;
	    void *  _value;
	    int     _size;
	    friend class Param;
};

class Param {
    public:
	    string extra;   /* location of the extra tpl file */

	    /* Energy calculation */
	    float scale_vdw0;
	    float scale_vdw1;
	    float scale_vdw;
	    float scale_tor;
	    float scale_ele;
	    float scale_dsolv;

	    /* Monte Carlo */
	    bool  monte_adding_conf;
	    bool  monte_average_pairwise;
	    float monte_warn_pairwise;
	    int   monte_seed;   /* random number seed */
	    float monte_temp;   /* temperature */
	    int   monte_flips;  /* maximum number of residues to flip in each MC step */
	    int   monte_nstart;
	    int   monte_neq;
	    float monte_reduce;
	    int   monte_runs;
	    int   monte_niter;
	    int   monte_trace;
	    int   nstate_max;

	    int   monte_n_red;
	    int   monte_niter_cycle;
	    int   monte_niter_min;
	    int   monte_niter_max;
	    int   monte_niter_chk;
	    float monte_converge;
	    int   monte_print_nonzero;

	    TitraType   titra_type;
	    int     titra_steps;
	    float   titra_ph0;
	    float   titra_phd;
	    float   titra_eh0;
	    float   titra_ehd;

	    float big_pairwise;

	    Param();
	    Param& load_env(const char*);
	    Param& load_env(ifstream&);
	    Param& load_extra(const char*);
	    Param& load_extra(ifstream&);
	    /* for the parameter database */
	    Param& load_param(ifstream& param_file);
	    void save(string key1, string key2, string key3, void *value, int size);
	    bool exist(string, string, string);
	    bool get(string key1, string key2, string key3);

    private:
	    map <string, Param_value> _db;
};

class Atom {
    friend class Prot;
    friend class Res;
    friend class Conf;

    public:
	    Prot*   prot;
	    Conf*   conf;
	    Math3D::Vector4<float> r;  /* atom position */
	    float   rad;
	    float   crg;
	    float   occ;
	    float   bFactor;
	    float   vdw_rad;
	    float   vdw_eps;
	    string  atomName;
	    string  atomType;

	    Atom&   read_pdbstr(const string&);
	    Atom&   read_orig_pdbstr(const string&);
	    string& write_pdbstr(string&);
	    Atom&   read_demetristr(const string&);
	    Atom&   switch_on()  {_on = true; return *this;};
	    Atom&   switch_off() {_on = false; return *this;};
	    bool    on()  {return _on;};
	    bool    off() {return !_on;};
	    bool    inRes(const Res&);
	    bool    inConf(const Conf&);
	    bool    isbkb();

    private:
	    /* Atom */
	    bool    _on;
	    int     atomSerial;
	    HetType  hetType;

	    /* conformer */
	    string  confName;
	    string  confID;
	    char    altLoc;
	    string  confHistory; /* 2 char for confName, Rotate, Swing, OH */

	    /* residue */
	    string  resName;
	    char    chainID;
	    int     resSeq;
	    char    iCode;

	    float   sas;

	    int    i_atom_prot;
	    int    i_atom_conf;
	    int    i_conf_res;
	    int    i_res_prot;
	    int    i_elem;
};

class Conf {
    friend class Prot;
    friend class Res;
    friend class Atom;
    friend bool sameID(Conf &conf1, Conf &conf2);

    public:
	    bool    on;
	    vector <Atom>  atom;
	    Prot*   prot;
	    Res*    res;
	    string  uniqID;
	    vector<float>  occ_table;
	    vector<float>  TS_table;
	    int    i_sidechain_prot;  /* global sidechain index in the protein */

	    Conf();
	    Conf(const Res&);
	    Conf(const Conf&);
	    Conf(const Atom&);

	    bool    inRes(const Res& _res);
	    int     ins_atom(const Atom&);
	    Conf&   initialize_by_uniqID(const string& _uniqID);
	    int     find_katom(const string& _atomName);
	    string  mk_uniqID();

	    /* residue */
	    string  resName;
	    char    chainID;
	    int     resSeq;
	    char    iCode;

	    /* conformer */
	    int     iConf;
	    string  confName;  /* conforname is resName+1st+2nd char in "history" */
	    string  confID;
	    char    altLoc;

	    float  netcrg;       /* net charge */
	    float  pKa;          /* pKa of the acid: HA <=> H+ + A- */
	    float  Em;           /* Em of the reduced: RED <=> OX + e- */
	    float  E_torsion;    /* torsion energy */
	    float  E_vdw0;       /* VDW potential to its own (including backbone atoms) */
	    float  E_vdw1;       /* VDW potential to backbone atoms (excluding its own) */
	    float  E_epol;       /* elec to backbone */
	    float  E_tors;
	    float  E_rxn0;
	    float  E_rxn;
	    float  E_extra;      /* extra energy to make the calc. match training exep. */
	    float  E_dsolv;
	    float  E_ph;
	    float  E_eh;
	    float  E_mfe;
	    float  E_base;
	    float  E_self0;      /* self energy without mfe */
	    double E_self;       /* total self energy, with mfe */
	    double TS;           /* entropy correction */
	    float  occ;          /* occupance */
	    float  occ_old;      /* occupance backup */
	    float  e;            /* number of electron(s) gained */
	    float  H;            /* number of proton(s) gained */
	    string   confHistory;  /* history of making this conformer */
	    int    counter;      /* General counter */
	    int    n_atom;
	    char   toggle;

	    /* back index */
	    int    i_res_prot;
	    int    i_conf_res;
	    int    i_conf_subres;

	    /* Monte Carlo */
	    int     i_sidechain_prot_orig;
	    int     counter_trial;
	    int     counter_accept;

    private:
	    int    subresID;
	    int    i_subres_res;
};

class Res {
    friend class Prot;
    friend class Conf;
    friend class Atom;

    public:
	    vector <Conf> conf;
	    Prot*   prot;
	    int    i_res_prot;
	    Res();
	    Res(const string&);
	    Res(const Atom&);
	    Res(const Conf&);
	    Conf& assignbkb();
	    int ins_conf(const Conf&);
	    int ins_conf(const Conf&, const int);
	    int del_conf(const int&);
	    string resID();
	    /*string resID(); -> err: undefined reference to `Res::resID()' @@@@@@@@@@@@@@@@@@@@@*/
	    void make_sumcrg();
	    string   resName;   /* used for identifying a residue */
	    char    chainID;      /* used for identifying a residue */
	    int     resSeq;       /* used for identifying a residue */
	    float   bFactor;

    /* For Monte Carlo */
    vector <Res*> ngh;
    float  fixed_occ;
    int     n_flip_max; /* Maximum number of residues to flip in one MC step */
    int     counter_trial;  /* Counter for number of MC trials made on this residue */

    Conf*   conf_w;     /* Pointer to the active conformer in the current microstate */
    Conf*   conf_new;   /* Pointer to the conformer MC trys to flip to */
    Conf*   conf_old;   /* Pointer to the conformer MC flips from */

    bool    print_sumcrg;
    vector<float>  sum_H;
    vector<float>  sum_e;
    vector<float>  sum_crg;

    vector< vector<int> > conf_type;

    private:
	    bool    on;
	    HetType hetType;
	    char    iCode;        /* used for identifying a residue */
	    int     cal_vdw;
	    bool    do_rot;       /* flag to do rotamers */
	    bool    do_sw;
	    bool    opt_hyd;
	    int     rotations;
	    float   phi_swing;
	    int     relax;
	    int     n_hyd_pos;      /* number of position for each degree of freedom */
	    float   sas;

	    float r12sq_max;
	    float r13sq_max;
	    float r14sq_max;

	    int    i_conf_on;
	    int    nconf_limit;
	    int    n_conf_ori;

    friend ostream& operator << (ostream& _file, Res& _res);
};

class Prot {
    friend class Res;
    friend class Conf;
    friend class Atom;

    public:
	    /* environment vairables */
	    Param   param;
	    vector <Res> res;
	    vector <Conf*> sidechain;

	    /* Monte Carlo */
	    vector <int> weighted_res_list;
	    double E_base;
	    double E_state;
	    double E_min;
	    double E_accum;
	    double E_free_unfold;

	    /* methods */
	    int ins_res(const Res&);
	    int del_res(const int&);
	    void index_off() {_index_flag = false;};
	    void index_on();
	    bool index_valid() { return _index_flag; };
	    int  count_conf();
	    int  count_atom();
	    Prot&  load_pdb(ifstream&);
	    Prot&  load_pdb(const char*);
	    Prot&  load_orig_pdb(ifstream&);
	    Prot&  load_orig_pdb(const char*);
	    Prot&  load_demetri(ifstream&);
	    Prot&  load_demetri(const char*);
	    Prot&  load_conflist(ifstream&);
	    Prot&  load_conflist(const char*);
	    Prot&  load_occtable(ifstream&);
	    Prot&  load_occtable(const char*);
	    Prot&  all_info2atom();
	    Prot&  write_pdb(ofstream&);
	    Prot&  write_pdb(ostream& pdb_file);
	    Prot&  write_pdb(const char*);
	    Prot&  write_conflist(ofstream&);
	    Prot&  write_conflist(const char*);
	    Res*   find_res(const Res&);
	    void   make_sum_crg();

	    vector <float>  titra_point;
	    vector <float>  sum_H;
	    vector <float>  sum_e;
	    vector <float>  sum_crg;

	    /* this should be temporary */
	    vector < vector<double> > _pairwise;

    private:
	    bool    _index_flag;
	    vector < vector<double> > _pairwise_lj;
	    vector < vector<double> > _pairwise_elec;
	    vector < vector<double> > _pairwise_active;
};

class Pka {
    public:
    bool    valid;
    bool    out_of_range;
    float   pK;
    float   n;
    float   chi2;
    float   coef; /* fitting coefficient */
};


class Pairwise {
    float ele;
    float vdw;
    float crt;
    float ori;
    char  mark[4];
};


typedef struct {
    char   uniqID[15];
    char   on;
    float  occ;          /* occupance */
    float  netcrg;       /* net charge */
    float  Em;           /* Em of the reduced: RED <=> OX + e- */
    float  pKa;          /* pKa of the acid: HA <=> H+ + A- */
    int    e;            /* number of electron(s) gained */
    int    H;            /* number of proton(s) gained */
    float  E_vdw0;       /* VDW potential to its own (including backbone atoms) */
    float  E_vdw1;       /* VDW potential to backbone atoms (excluding its own) */
    float  E_tors;
    float  E_epol;       /* elec to backbone */
    float  E_rxn0;
    float  E_rxn;
    float  E_dsolv;
    float  E_extra;      /* extra energy to make the calc. match training exep. */
    char   history[10];  /* history of making this conformer */
} CONF_HEAD;

class Ematrix {
    vector <Conf> conf;
    Pairwise  **pw;   
};

typedef Math3D::Vector4<float>  Vect;
typedef vector <Atom>::iterator Atom_iter;
typedef vector <Conf>::iterator Conf_iter;
typedef vector  <Res>::iterator Res_iter;

void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk);
Pka curve_fitting(vector <float> titra_point, vector <float> sumcrg_table, float shift_pK, float shift_n);
Pka curve_fitting(vector <float> titra_point, vector <float> sumcrg_table);
vector<Pka> curve_fitting_2pK(vector <float> titra_point, vector <float> sumcrg_table);
vector<Pka> curve_fitting_N_pK(vector <float> titra_point, vector <float> sumcrg_table, int& N_pK);

/* extra string functions */
string rm_comment(const string& str);
string strip(const string& str);
string strip_tail(const string& str);
vector<string> split(const string& str);
vector<string> split(const string& str, const char* delimiters);
vector<string> file2str(ifstream& infile);
vector<string> file2str(const char* file_name);
bool equ_str(const string& str1, const string& str2);
vector <string> get_files(const char *dir_name);

int load_opp(vector <double>&, ifstream&, const Prot&, const int);
int load_opp(vector <double>&, const char*, const Prot&, const int);
int load_opp(double*, ifstream&, const Prot&);
int load_opp(double*, const char*, const Prot&);

/* extra vector functions */
vector<float> sum_2vec(vector<float> vec1, vector<float> vec2);
void print_vec(const vector<float> vec, const string fmt);
string vec2string(const vector<float> vec, const string fmt);
string vec2string(const vector<int> vec, const string fmt);
vector <float> get_stats(vector<float> list, vector<float> weight, float bin_size);
vector <float> get_stats(vector<float> list, float bin_size);
vector < vector <int> > get_bins(vector<int> list, vector<float> values, float bin_size);
vector < vector <int> > get_bins(vector<int> list, vector<float> values, float bin_size, float zero_point);
float get_dihedral(Vector4<float> v0, Vector4<float> v1, Vector4<float> v2, Vector4<float> v3);


ostream& operator << (ostream& _file, Res& _res);

/* methods for a vector */
template<class T> T sum(const vector<T>& values)
{
    T sum = 0.;
    for (int i=0;i<values.size();++i) sum += values[i];
    return sum;
}

template<class T> T sumsq(const vector<T>& values)
{
    T sumsq = 0.;
    for (int i=0;i<values.size();++i) sumsq += values[i]*values[i];
    return sumsq;
}

template<class T> T average(const vector<T>& values)
{
    if (values.size())
        return sum(values)/((T) values.size());
    else
        return 0.;
}

template<class T> T stdev(const vector<T>& values)
{
    if (values.size() <= 1) return 0.0;
    else {
        T tmp  = sum(values);
        T tmp2 = (sumsq(values) - tmp*tmp/((double)values.size()))/((double)(values.size()-1));
        if (tmp2 <= 0.) return 0.0; /* avoid precision error */
        else return sqrt(tmp2);
    }
}

template<class T> T min(const vector<T>& values)
{
    T min_val;
    int i;
    if (values.size()) {
        min_val = values[0];
    }
    else {
        min_val = 0;
    }
    
    for (i=1; i<values.size(); ++i) {
        if (min_val > values[i])
            min_val = values[i];
    }
    
    return min_val;
}

template<class T> T max(const vector<T>& values)
{
    T max_val;
    int i;
    if (values.size()) {
        max_val = values[0];
    }
    else {
        max_val = 0;
    }
    
    for (i=1; i<values.size(); ++i) {
        if (max_val < values[i])
            max_val = values[i];
    }
    
    return max_val;
}

template<class T> T max_abs(const vector<T>& values)
{
    T max_val = 0;
    for (int i=0; i<values.size(); ++i) {
        if (max_val < fabs(values[i]))
            max_val = fabs(values[i]);
    }
    return max_val;
}

template<class T> int min_index(const vector<T>& values)
{
    T min_val;
    int i, i_min;
    if (values.size()) {
        min_val = values[0];
    }
    else {
        min_val = 0;
    }
    
    for (i=1; i<values.size(); ++i) {
        if (min_val > values[i]) {
            min_val = values[i];
            i_min = i;
        }
    }
    
    return i_min;
}

template<class T> int nearest_index(const T search_val, const vector<T>& values)
{
    if (!values.size()) return -1;
    T min_diff = fabs(search_val - values[0]);
    int ret_idx = 0;
    
    for (int i=0; i<values.size(); ++i) {
        if ( fabs(search_val - values[i]) < min_diff ) {
            min_diff = fabs(search_val - values[i]);
            ret_idx = i;
        }
    }
    return ret_idx;
}

template<class T> int find_match_index(const T search_val, const vector<T>& values, const T threshold)
{
    for (int i=0; i<values.size(); ++i) {
        if ( fabs(search_val - values[i]) < threshold ) {
            return i;
        }
    }
    return -1;
}

template<class T> int rfind_match_index(const T search_val, const vector<T>& values, const T threshold)
{
    for (int i=values.size()-1; i>=0; --i) {
        if ( fabs(search_val - values[i]) < threshold ) {
            return i;
        }
    }
    return -1;
}

template<class T> vector<T> subvec(const vector<T>& vec, const int idx_start, const int idx_end)
{
    vector<T> ret_val;
    int _idx_start = idx_start;
    int _idx_end   = idx_end;
    if (_idx_start < 0) _idx_start = 0;
    if (_idx_end >= vec.size()) _idx_end = vec.size()-1;
    for (int i=_idx_start; i<=_idx_end; ++i) {
        ret_val.push_back(i);
    }
    return ret_val;
}

/* Numerical Recipes */
double ran2(long*);

/* Modules */
int Monte();
