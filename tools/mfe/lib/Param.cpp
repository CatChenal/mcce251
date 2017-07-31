#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

Param::Param()
{
    scale_ele     = 1.0;
    scale_dsolv   = 1.0;
    scale_tor     = 1.0;
    scale_vdw0    = 1.0;
    scale_vdw1    = 1.0;
    scale_vdw     = 1.0;
}

Param& Param::load_env(const char* file_name)
{
    ifstream infile;
    string msg;

    infile.open(file_name, ios::in);
    if (!infile) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_env(infile);
    infile.close();
    return *this;
}

Param& Param::load_env(ifstream& infile)
{
    string line;

    /* initialize */
    scale_ele     = 1.0;
    scale_dsolv   = 1.0;
    scale_tor     = 1.0;
    scale_vdw0    = 1.0;
    scale_vdw1    = 1.0;
    scale_vdw     = 1.0;

    monte_average_pairwise = true;
    monte_warn_pairwise = 200.;
    monte_niter_max   =   -1;
    monte_niter_chk   =  100;
    monte_converge    = 1e-4;
    monte_print_nonzero =  1;

    monte_seed   = -1;
    monte_temp   = 298.15;
    monte_nstart = 100;
    monte_n_red  = -1;
    monte_neq    = 200;
    monte_niter  = 5000;
    monte_flips  = 3;
    monte_reduce = 0.001;
    monte_niter_cycle = 500000;
    monte_niter_min = 1000;
    monte_niter_max = 5000;
    monte_niter_chk = 500;
    monte_converge = 0.05;

    while(getline(infile, line)) {
        line = strip(rm_comment(line));

        if (line.find("(EXTRA)")!=string::npos) {
            string str1 = line.substr(0,line.find(' ')).c_str();
            Param::extra = str1;
        }
        /* Monte Carlo */
        else if (line.find("(MONTE_SEED)")!=string::npos) {
            Param::monte_seed = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_T)")!=string::npos) {
            Param::monte_temp = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_RUNS)")!=string::npos) {
            Param::monte_runs = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NSTART)")!=string::npos) {
            Param::monte_nstart = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NEQ)")!=string::npos) {
            Param::monte_neq = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NITER)")!=string::npos) {
            Param::monte_niter = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_FLIPS)")!=string::npos) {
            Param::monte_flips = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_TRACE)")!=string::npos) {
            Param::monte_trace = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_N_REDUCE)")!=string::npos) {
            Param::monte_n_red = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_REDUCE)")!=string::npos) {
            Param::monte_reduce = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NITER_CYCLE)")!=string::npos) {
            Param::monte_niter_cycle = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NITER_MIN)")!=string::npos) {
            Param::monte_niter_min = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NITER_MAX)")!=string::npos) {
            Param::monte_niter_max = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_NITER_CHK)")!=string::npos) {
            Param::monte_niter_chk = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_CONVERGE)")!=string::npos) {
            Param::monte_converge = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(MONTE_PRINT_NONZERO)")!=string::npos) {
            string str1 = line.substr(0,line.find(' ')).c_str();
            if (str1[0] == 't' || str1[0] == 'T') {
                Param::monte_print_nonzero = 1;
            }
            else Param::monte_print_nonzero = 0;
        }
        else if (line.find("(TITR_TYPE)")!=string::npos) {
            string str1 = line.substr(0,line.find(' ')).c_str();
            if (str1[0] == 'p' || str1[0] == 'P') Param::titra_type = PH;
            else Param::titra_type = EH;
        }
        else if (line.find("(TITR_PH0)")!=string::npos) {
            Param::titra_ph0 = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(TITR_EH0)")!=string::npos) {
            Param::titra_eh0 = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(TITR_PHD)")!=string::npos) {
            Param::titra_phd = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(TITR_EHD)")!=string::npos) {
            Param::titra_ehd = atof(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(TITR_STEPS)")!=string::npos) {
            Param::titra_steps = atoi(line.substr(0,line.find(' ')).c_str());
        }
        else if (line.find("(BIG_PAIRWISE)")!=string::npos) {
            Param::big_pairwise = atof(line.substr(0,line.find(' ')).c_str());
        }
    }
    return *this;
}

Param& Param::load_extra(const char* file_name)
{
    ifstream infile;
    string msg;

    infile.open(file_name, ios::in);
    if (!infile) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_extra(infile);
    infile.close();
    return *this;
}

Param& Param::load_extra(ifstream& infile)
{
    string line;
    while(getline(infile, line)) {
        if (line.find("SCALING ")!=string::npos) {
            if (line.find(" ELE ")!=string::npos) {
                scale_ele = atof(line.substr(15).c_str());
            }
            else if (line.find(" DSOLV ")!=string::npos) {
                scale_dsolv = atof(line.substr(15).c_str());
            }
            else if (line.find(" VDW0 ")!=string::npos) {
                scale_vdw0 = atof(line.substr(15).c_str());
            }
            else if (line.find(" VDW1 ")!=string::npos) {
                scale_vdw1 = atof(line.substr(15).c_str());
            }
            else if (line.find(" VDW ")!=string::npos) {
                scale_vdw = atof(line.substr(15).c_str());
            }
            else if (line.find(" TORS ")!=string::npos) {
                scale_tor = atof(line.substr(15).c_str());
            }
        }
    }
    return *this;
}

void Param::save(string key1, string key2, string key3, void *value, int size) {
    string  key;
    Param_value param;

    key = strip(key1) + strip(key2) + strip(key3);

    param._key = strip(key1) + strip(key2) + strip(key3);
    param._value = value;
    param._size = size;

    _db[key] = param;
}

bool Param::exist(string key1, string key2, string key3) {
    string  key;
    key = strip(key1) + strip(key2) + strip(key3);

    if (_db.find(key) == _db.end()) {
        return false;
    }
    else {
        return true;
    }
}

bool Param::get(string key1, string key2, string key3) {
    string  key;
    key = strip(key1) + strip(key2) + strip(key3);
    if (!exist(key1,key2,key3)) return false;

    if (_db.find(key) == _db.end()) {
        return false;
    }
    else {
        return true;
    }
}

#define LEN_KEY1 10
#define LEN_KEY2 5
#define LEN_KEY3 5

typedef struct {
    bool ligand;
    int  res_offset;
    string name;
} Connected_atom;

typedef struct {
    string orbital;
    vector <Connected_atom> atom;
} Connectivity;

int load_param(ifstream& param_file)
{
    string line;

    while(getline(param_file, line)) {
        line = strip_tail(rm_comment(line));
        if (line.size() < 20) continue;

        string key1 = strip(line.substr(0,LEN_KEY1));
        string key2 = strip(line.substr(LEN_KEY1,LEN_KEY2));
        string key3 = strip(line.substr(LEN_KEY1+LEN_KEY2,LEN_KEY3));

        if (key1 == "CONNECT") {
            Connectivity* connect = new Connectivity;
            connect->orbital = strip(line.substr(LEN_KEY1+LEN_KEY2+LEN_KEY3, 9));
            int i_substr = LEN_KEY1+LEN_KEY2+LEN_KEY3+10;
            while (i_substr < line.size()) {
                Connected_atom connected;
                connected.ligand = ( strip(line.substr(i_substr,4)) == "LIG" ) ? true:false;
                connected.res_offset = atoi(line.substr(i_substr,4).c_str());
                connected.name = line.substr(i_substr+5,4);
                while (connected.name.size() < 4) connected.name += " ";
                connect->atom.push_back(connected);
                i_substr += 10;
            }
        }
        else {
        }
    }
    return 0;
}
