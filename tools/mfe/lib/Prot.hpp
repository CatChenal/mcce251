#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

inline int Prot::ins_res(const Res& _res)
{
    index_off();
    res.push_back(_res);
    res.back().prot = this;
    
    return res.size()-1;
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
            res[k_res].asignbkb();
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
    
    return *this;
}

Prot& Prot::load_pdb(char* file_name)
{
    ifstream inpdb;
    string msg;
    
    inpdb.open(file_name, ios::in);
    if (!inpdb) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    load_pdb(inpdb);
    return *this;
}

Prot& Prot::write_pdb(ofstream& pdb_file)
{
    string pdbstr;
    Res_iter i_res;
    Conf_iter i_conf;
    Atom_iter i_atom;
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

Prot& Prot::write_pdb(char* file_name)
{
    ofstream outpdb;
    string msg;
    
    outpdb.open(file_name, ios::out);
    if (!outpdb) {
        cerr << "Error: can't open file " << file_name << " for writing." << endl;
        throw msg = "write_file_error";
    }
    write_pdb(outpdb);
    return *this;
}

Prot& Prot::load_conflist(ifstream& conflist)
{
    string line;
    vector<string> tokens;
    int i_column;
    
    /* skip first line */
    getline(conflist, line);
    
    while(getline(conflist, line)) {
        int k_res, k_conf;
        Conf temp_conf;
        if (line.size() < 20) continue;
        
        string uniqID = line.substr(6,14);
        temp_conf.initialize_by_uniqID(uniqID);
        temp_conf.toggle    = line[21];
        temp_conf.occ       = atof(line.substr(23, 4).c_str());
        temp_conf.netcrg    = atof(line.substr(28, 6).c_str());
        temp_conf.Em        = atof(line.substr(35, 5).c_str());
        temp_conf.pKa       = atof(line.substr(41, 5).c_str());
        temp_conf.e         = atoi(line.substr(47, 2).c_str());
        temp_conf.H         = atoi(line.substr(50, 2).c_str());
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
            res[k_res].asignbkb();
        }
        
        k_conf = res[k_res].ins_conf(temp_conf);
    }
    
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
        titraType = ph;
    }
    else if (tokens[0] == "eh") {
        titraType = eh;
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
            res[k_res].asignbkb();
        }
        
        /* find the matching residue */
        for (k_conf = res[k_res].conf.size() -1; k_conf >= 1; --k_conf) {
            if (temp_conf.confName == res[k_res].conf[k_conf].confName &&
                temp_conf.confID   == res[k_res].conf[k_conf].confID) break;
        }
        if (k_conf < 1) {
            k_conf = res[k_res].ins_conf(temp_conf);
        }
        
        for (i_column = 1; i_column < tokens.size(); ++i_column) {
            res[k_res].conf[k_conf].occ_table.push_back( atof( tokens[i_column].c_str() ) );
        }
    }
    
    return *this;
}

inline int Prot::count_atom()
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

inline void Prot::index_on() {
    if (index_valid()) return;
    int i_res, i_conf, i_atom;
    
    sidechain.clear();
    for (i_res = 0; i_res < res.size(); ++i_res) {
        for (i_conf = 1; i_conf < res[i_res].conf.size(); ++i_conf) {
            res[i_res].conf[i_conf].i_conf_prot = sidechain.size();
            sidechain.push_back(&res[i_res].conf[i_conf]);
        }
    }
    _index_flag = true;
    return;
}
