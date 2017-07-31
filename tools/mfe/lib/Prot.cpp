#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

int Prot::ins_res(const Res& _res)
{
    index_off();
    res.push_back(_res);
    res.back().prot = this;

    return res.size()-1;
}

int Prot::del_res(const int& i_res)
{
    if (i_res<0) return -1;
    if (i_res >= res.size()) return -1;

    index_off();
    res.erase(res.begin()+i_res);
    return 0;
}

Prot& Prot::load_pdb(ifstream& pdb_file)
{
    string pdbline;

    while(getline(pdb_file, pdbline)) {
        Atom      atom;
        int       k_res, k_conf;

        atom.read_pdbstr(pdbline);
        if (atom.off()) continue;

        /* find the matching residue */
        for (k_res = res.size() -1; k_res >= 0; --k_res) {
            if (atom.inRes(res[k_res])) break;
        }
        /* If this atom doesn't belong to any existing residue, add a new residue */
        if (k_res < 0) {
            Res temp_res(atom);
            k_res = ins_res(temp_res);
            res[k_res].assignbkb();
        }

        if (atom.isbkb()) k_conf = 0;
        else {
            for (k_conf = res[k_res].conf.size()-1; k_conf >= 0; --k_conf) {
                if (atom.inConf(res[k_res].conf[k_conf])) break;
            }
            /* If this atom doesn't belong to any existing conformer, add a new conformer */
            if (k_conf < 0) {
                Conf temp_conf(atom);
                k_conf = res[k_res].ins_conf(temp_conf);
            }
        }

        /* Add atom */
        res[k_res].conf[k_conf].ins_atom(atom);
    }

    for (int i_res=0; i_res<res.size(); ++i_res) {
        for (int i_conf=0; i_conf<res[i_res].conf.size(); ++i_conf) {
            res[i_res].conf[i_conf].mk_uniqID();
        }
    }
    return *this;
}

Prot& Prot::load_pdb(const char* file_name)
{
    ifstream inpdb;
    string msg;

    inpdb.open(file_name, ios::in);
    if (!inpdb) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        return *this;
        //throw msg = "read_file_error";
    }
    load_pdb(inpdb);
    return *this;
}

Prot& Prot::load_orig_pdb(ifstream& pdb_file)
{
    string pdbline;

    while(getline(pdb_file, pdbline)) {
        Atom      atom;
        int       k_res, k_conf;

        atom.read_orig_pdbstr(pdbline);
        if (atom.off()) continue;

        /* find the matching residue */
        for (k_res = res.size() -1; k_res >= 0; --k_res) {
            if (atom.inRes(res[k_res])) break;
        }
        /* If this atom doesn't belong to any existing residue, add a new residue */
        if (k_res < 0) {
            Res temp_res(atom);
            k_res = ins_res(temp_res);
            res[k_res].assignbkb();
        }

        if (atom.isbkb()) k_conf = 0;
        else {
            for (k_conf = res[k_res].conf.size()-1; k_conf >= 0; --k_conf) {
                if (atom.inConf(res[k_res].conf[k_conf])) break;
            }
            /* If this atom doesn't belong to any existing conformer, add a new conformer */
            if (k_conf < 0) {
                Conf temp_conf(atom);
                k_conf = res[k_res].ins_conf(temp_conf);
            }
        }
        /* Add atom */
        res[k_res].conf[k_conf].ins_atom(atom);
    }

    for (int i_res=0; i_res<res.size(); ++i_res) {
        for (int i_conf=0; i_conf<res[i_res].conf.size(); ++i_conf) {
            res[i_res].conf[i_conf].mk_uniqID();
        }
    }
    return *this;
}

Prot& Prot::load_orig_pdb(const char* file_name)
{
    ifstream inpdb;
    string msg;

    inpdb.open(file_name, ios::in);
    if (!inpdb) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_orig_pdb(inpdb);
    return *this;
}

Prot& Prot::load_demetri(ifstream& pdb_file)
{
    string pdbline;

    while(getline(pdb_file, pdbline)) {
        Atom      atom;
        int       k_res, k_conf;

        atom.read_demetristr(pdbline);
        if (atom.off()) continue;

        /* find the matching residue */
        for (k_res = res.size() -1; k_res >= 0; --k_res) {
            if (atom.inRes(res[k_res])) break;
        }
        /* If this atom doesn't belong to any existing residue, add a new residue */
        if (k_res < 0) {
            Res temp_res(atom);
            k_res = ins_res(temp_res);
            res[k_res].assignbkb();
        }

        if (atom.isbkb()) k_conf = 0;
        else {
            for (k_conf = res[k_res].conf.size()-1; k_conf >= 0; --k_conf) {
                if (atom.inConf(res[k_res].conf[k_conf])) break;
            }
            /* If this atom doesn't belong to any existing conformer, add a new conformer */
            if (k_conf < 0) {
                Conf temp_conf(atom);
                k_conf = res[k_res].ins_conf(temp_conf);
            }
        }

        /* Add atom */
        res[k_res].conf[k_conf].ins_atom(atom);
    }

    for (int i_res=0; i_res<res.size(); ++i_res) {
        for (int i_conf=0; i_conf<res[i_res].conf.size(); ++i_conf) {
            res[i_res].conf[i_conf].mk_uniqID();
        }
    }
    return *this;
}

Prot& Prot::load_demetri(const char* file_name)
{
    ifstream inpdb;
    string msg;

    inpdb.open(file_name, ios::in);
    if (!inpdb) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_demetri(inpdb);
    return *this;
}

Prot& Prot::load_conflist(ifstream& conflist)
{
    string line;
    vector<string> tokens;

    /* skip first line */
    getline(conflist, line);

    while(getline(conflist, line)) {
        int k_res, k_conf;
        Conf temp_conf;
        line = rm_comment(line);
        if (line.size() < 20) continue;

        /*
        .012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        .iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    self
        .02082 _NA+1Y0010_002 f 0.00  1.000     0  0.00  0  0   0.000  -0.002   0.000  -0.023   0.678   4.000 +1O000M000t
        .00006 ARG+1A0014_003 f 0.00  1.000     0 12.50  0  0   1.000   1.318   6.424   3.089   0.639   0.232  -1.750
        */
        string uniqID = line.substr(6,14);
        temp_conf.initialize_by_uniqID(uniqID);
        temp_conf.iConf     = atoi(line.substr( 0, 5).c_str());
        temp_conf.toggle    = line[21];
        temp_conf.occ       = atof(line.substr(23, 4).c_str());
        temp_conf.netcrg    = atof(line.substr(28, 6).c_str());
        temp_conf.Em        = atof(line.substr(35, 5).c_str());
        temp_conf.pKa       = atof(line.substr(41, 5).c_str());
        temp_conf.e         = atof(line.substr(47, 2).c_str());
        temp_conf.H         = atof(line.substr(50, 2).c_str());
        temp_conf.E_vdw0    = atof(line.substr(53, 7).c_str());
        temp_conf.E_vdw1    = atof(line.substr(61, 7).c_str());
        temp_conf.E_tors    = atof(line.substr(69, 7).c_str());
        temp_conf.E_epol    = atof(line.substr(77, 7).c_str());
        temp_conf.E_dsolv   = atof(line.substr(85, 7).c_str());
        temp_conf.E_extra   = atof(line.substr(93, 7).c_str());
        temp_conf.E_self    = atof(line.substr(101,7).c_str());

        /* find the matching residue */
        for (k_res = res.size() -1; k_res >= 0; --k_res) {
            if (temp_conf.inRes(res[k_res])) break;
        }
        /* If couldn't find the residue, add a new one */
        if (k_res < 0) {
            Res temp_res(temp_conf);
            k_res = ins_res(temp_res);
            res[k_res].assignbkb();
        }

        k_conf = res[k_res].ins_conf(temp_conf);
    }
    return *this;
}

Prot& Prot::load_conflist(const char* file_name)
{
    ifstream conflist;
    string msg;

    conflist.open(file_name, ios::in);
    if (!conflist) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_conflist(conflist);
    return *this;
}

Prot& Prot::load_occtable(ifstream& occtable)
{
    string line;
    vector<string> tokens;
    int i_column;

    /* use first line to get ph/eh */
    getline(occtable, line);
    tokens = split(line);
    if (tokens[0] == "ph") {
        param.titra_type = PH;
    }
    else if (tokens[0] == "eh") {
        param.titra_type = EH;
    }
    for (i_column = 1; i_column < tokens.size(); ++i_column) {
        titra_point.push_back( atof( tokens[i_column].c_str() ) );
    }

    while(getline(occtable, line)) {
        int k_res, k_conf;
        Conf temp_conf;
        if (line.size() < 14) continue;

        tokens = split(line);

        string uniqID = line.substr(0,14);
        temp_conf.initialize_by_uniqID(uniqID);

        /* find the matching residue */
        for (k_res = res.size() -1; k_res >= 0; --k_res) {
            if (temp_conf.inRes(res[k_res])) break;
        }
        /* If couldn't find the residue, add a new one */
        if (k_res < 0) {
            Res temp_res(temp_conf);
            k_res = ins_res(temp_res);
            res[k_res].assignbkb();
        }

        /* find the matching conformer */
        for (k_conf = res[k_res].conf.size() -1; k_conf >= 1; --k_conf) {
            if (temp_conf.confName == res[k_res].conf[k_conf].confName &&
                temp_conf.confID   == res[k_res].conf[k_conf].confID) break;
        }
        /* If couldn't find the conformer, add a new one */
        if (k_conf < 1) {
            k_conf = res[k_res].ins_conf(temp_conf);
        }

        for (i_column = 1; i_column < tokens.size(); ++i_column) {
            res[k_res].conf[k_conf].occ_table.push_back( atof( tokens[i_column].c_str() ) );
        }
    }
    return *this;
}

Prot& Prot::load_occtable(const char* file_name)
{
    ifstream occtable;
    string msg;
    occtable.open(file_name, ios::in);
    if (!occtable) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_occtable(occtable);
    occtable.close();
    return *this;
}

Prot& Prot::all_info2atom()
{
    for (int i_res = 0; i_res < res.size(); ++i_res) {
        for (int i_conf = 0; i_conf < res[i_res].conf.size(); ++i_conf) {
            for (int i_atom = 0; i_atom < res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                res[i_res].conf[i_conf].atom[i_atom].hetType = res[i_res].hetType;
                res[i_res].conf[i_conf].atom[i_atom].chainID = res[i_res].chainID;
                res[i_res].conf[i_conf].atom[i_atom].resName = res[i_res].resName;
                res[i_res].conf[i_conf].atom[i_atom].resSeq  = res[i_res].resSeq;
                res[i_res].conf[i_conf].atom[i_atom].iCode   = res[i_res].iCode;

                res[i_res].conf[i_conf].atom[i_atom].confName    = res[i_res].conf[i_conf].confName;
                res[i_res].conf[i_conf].atom[i_atom].confID      = res[i_res].conf[i_conf].confID;
                res[i_res].conf[i_conf].atom[i_atom].altLoc      = res[i_res].conf[i_conf].altLoc;
                res[i_res].conf[i_conf].atom[i_atom].confHistory = res[i_res].conf[i_conf].confHistory;
            }
        }
    }
    return *this;
}

Prot& Prot::write_pdb(ofstream& pdb_file)
{
    string pdbstr;
    Res_iter i_res;
    Conf_iter i_conf;
    Atom_iter i_atom;
    all_info2atom();
    for (i_res = res.begin(); i_res < res.end(); ++i_res) {
        for (i_conf = i_res->conf.begin(); i_conf < i_res->conf.end(); ++i_conf) {
            for (i_atom = i_conf->atom.begin(); i_atom < i_conf->atom.end(); ++i_atom) {
                i_atom->write_pdbstr(pdbstr);
                pdb_file << pdbstr << endl;
            }
        }
    }
    return *this;
}

Prot& Prot::write_pdb(const char* file_name)
{
    ofstream outpdb;
    string msg;

    if (strlen(file_name) == 0) {
        write_pdb(cout);
    }
    else {
        outpdb.open(file_name, ios::out);
        if (!outpdb) {
            write_pdb(cout);
        }
        else {
            write_pdb(outpdb);
        }
    }
    return *this;
}

Prot& Prot::write_pdb(ostream& pdb_out)
{
    string pdbstr;
    Res_iter i_res;
    Conf_iter i_conf;
    Atom_iter i_atom;
    all_info2atom();
    for (i_res = res.begin(); i_res < res.end(); ++i_res) {
        for (i_conf = i_res->conf.begin(); i_conf < i_res->conf.end(); ++i_conf) {
            for (i_atom = i_conf->atom.begin(); i_atom < i_conf->atom.end(); ++i_atom) {
                i_atom->write_pdbstr(pdbstr);
                pdb_out << pdbstr << endl;
            }
        }
    }
    return *this;
}

Prot& Prot::write_conflist(ofstream& conflist)
{
    char cstring[160];

    conflist << "iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    self\n";
    for (int i_res = 0; i_res < res.size(); ++i_res) {
        for (int i_conf = 1; i_conf < res[i_res].conf.size(); ++i_conf) {

            sprintf(cstring, "%05d %s %c %4.2f %6.3f %5.0f %5.2f %2.0f %2.0f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f",
            (res[i_res].conf[i_conf].iConf <= 99999 ? res[i_res].conf[i_conf].iConf : 99999),
            res[i_res].conf[i_conf].uniqID.c_str(),
            res[i_res].conf[i_conf].toggle,
            res[i_res].conf[i_conf].occ,
            res[i_res].conf[i_conf].netcrg,
            res[i_res].conf[i_conf].Em,
            res[i_res].conf[i_conf].pKa,
            res[i_res].conf[i_conf].e,
            res[i_res].conf[i_conf].H,
            res[i_res].conf[i_conf].E_vdw0,
            res[i_res].conf[i_conf].E_vdw1,
            res[i_res].conf[i_conf].E_tors,
            res[i_res].conf[i_conf].E_epol,
            res[i_res].conf[i_conf].E_dsolv,
            res[i_res].conf[i_conf].E_extra,
            res[i_res].conf[i_conf].E_self );

            conflist << cstring << endl;
        }
    }
    return *this;
}

Prot& Prot::write_conflist(const char* file_name)
{
    ofstream conflist;
    string msg;

    conflist.open(file_name, ios::out);
    if (!conflist) {
        cerr << "Error: can't open file " << file_name << " for writing." << endl;
        throw msg = "write_file_error";
    }
    write_conflist(conflist);
    return *this;
}

int Prot::count_atom()
{
    int i_res, i_conf, i_atom;
    int counter = 0;
    for (i_res = 0; i_res < res.size(); ++i_res) {
        for (i_conf = 0; i_conf < res[i_res].conf.size(); ++i_conf) {
            for (i_atom = 0; i_atom < res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                if (res[i_res].conf[i_conf].atom[i_atom].off()) continue;
                counter++;
            }
        }
    }
    return counter;
}

int Prot::count_conf()
{
    int i_res, i_conf;
    int counter = 0;
    for (i_res = 0; i_res < res.size(); ++i_res) {
        for (i_conf = 0; i_conf < res[i_res].conf.size(); ++i_conf) {
            counter++;
        }
    }
    return counter;
}

void Prot::index_on()
{
    if (index_valid()) return;
    int i_res, i_conf, i_atom;

    /* i_sidechain_prot, sidechain */
    sidechain.clear();
    for (i_res = 0; i_res < res.size(); ++i_res) {
        res[i_res].prot = this;
        res[i_res].i_res_prot = i_res;
        for (i_conf = 0; i_conf < res[i_res].conf.size(); ++i_conf) {
            res[i_res].conf[i_conf].prot = this;
            res[i_res].conf[i_conf].res = &res[i_res];
            res[i_res].conf[i_conf].i_conf_res = i_conf;

            if (i_conf) {
                res[i_res].conf[i_conf].i_sidechain_prot = sidechain.size();
                sidechain.push_back(&res[i_res].conf[i_conf]);
            }

            for (i_atom = 0; i_atom < res[i_res].conf[i_conf].atom.size(); ++i_atom) {
                res[i_res].conf[i_conf].atom[i_atom].prot = this;
                res[i_res].conf[i_conf].atom[i_atom].conf = &res[i_res].conf[i_conf];
                res[i_res].conf[i_conf].atom[i_atom].i_atom_conf = i_atom;
            }
        }
    }
    _index_flag = true;
    return;
}

void Prot::make_sum_crg()
{
    sum_H.clear();
    sum_H.resize(titra_point.size(), 0);
    sum_e.clear();
    sum_e.resize(titra_point.size(), 0);
    sum_crg.clear();
    sum_crg.resize(titra_point.size(), 0);

    for (int i_res=0; i_res<res.size(); ++i_res) {
        res[i_res].make_sumcrg();

        sum_H   = sum_2vec(res[i_res].sum_H, sum_H);
        sum_e   = sum_2vec(res[i_res].sum_e, sum_e);
        sum_crg = sum_2vec(res[i_res].sum_crg, sum_crg);
    }
}

Res* Prot::find_res(const Res& search_res)
{
    int i_res;
    for (i_res = 0; i_res < res.size(); ++i_res) {
        char search_chainID     = search_res.chainID;
        char i_res_chianID      = res[i_res].chainID;
        char search_iCode       = search_res.iCode;
        char i_res_iCode        = res[i_res].iCode;

        if (search_chainID  == ' ') search_chainID  = '_';
        if (i_res_chianID   == ' ') i_res_chianID   = '_';
        if (search_iCode    == ' ') search_iCode    = '_';
        if (i_res_iCode     == ' ') i_res_iCode     = '_';

        if (   search_res.resName   == res[i_res].resName
            && search_chainID       == i_res_chianID
            && search_res.resSeq    == res[i_res].resSeq
            && search_iCode         == i_res_iCode ) break;
    }
    if (i_res >= res.size()) return NULL;
    else return &res[i_res];
}
