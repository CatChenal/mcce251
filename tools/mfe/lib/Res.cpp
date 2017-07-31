#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "mcce.hpp"

using namespace std;

bool sameID(Conf &conf1, Conf &conf2);

Res::Res()
{
};

Res::Res(const string& _resID)
{
    if (_resID.size() >= 9) {
        resName = _resID.substr(0,3);
        chainID = _resID[4];
        resSeq = atoi(_resID.substr(5,4).c_str());

        if (_resID.size() >= 10) {
            iCode = _resID[9];
        }
        else {
            iCode = '_';
        }
    }
}
Res::Res(const Atom& _atom)
{
    hetType    = _atom.hetType;
    chainID     = _atom.chainID;
    resName     = _atom.resName;
    resSeq      = _atom.resSeq;
    iCode       = _atom.iCode;
};

Res::Res(const Conf& _conf)
{
    chainID     = _conf.chainID;
    resName     = _conf.resName;
    resSeq      = _conf.resSeq;
    iCode       = _conf.iCode;
};

Conf& Res::assignbkb()
{
    Conf bkb(*this);
    bkb.confName =  resName+"BK";
    bkb.confID   = "000";
    bkb.confHistory  = "BK________";
    bkb.altLoc   = ' ';

    if (conf.empty()) {
        this->prot->index_off();
        conf.push_back(bkb);
        conf.back().prot = this->prot;
    }
    else {
        string err_msg;
        cerr << "   Error: Trying to add a backbone conformer into a non-empty conformer array.\n";
        throw err_msg = "conf_not_empty";
    }
    return conf[0];
}

bool sameID(Conf &conf1, Conf &conf2)
{
    if (conf1.confName == conf2.confName) {
        if (conf1.confID == conf2.confID) {
            if (conf1.resSeq == conf2.resSeq) {
                if (conf1.chainID == conf2.chainID) {
                    if (conf1.iCode == conf2.iCode) {
                        return true;
                    }
                    else
                        return false;
                }
                else
                    return false;
            }
            else
                return false;
        }
        else
            return false;
    }
    else
        return false;
}

int Res::ins_conf(const Conf& _conf)
{
    if (conf.empty()) {
        string err_msg;
        cerr << "   Error: This residue has no backbone yet.\n";
        throw err_msg = "no_backbone";
    }

    this->prot->index_off();
    conf.push_back(_conf);
    conf.back().prot = this->prot;

    return conf.size()-1;
}

int Res::ins_conf(const Conf& _conf, const int insert_pos)
{
    /* Error checking */
    if (conf.empty()) {
        string err_msg;
        cerr << "   Error: This residue has no backbone yet.\n";
        throw err_msg = "no_backbone";
    }

    if (insert_pos <= 0 || insert_pos > conf.size()) {
        string err_msg;
        cerr << "   Error: Cannot insert a conformer at position " << insert_pos << ".\n";
        throw err_msg = "out_of_bound";
    }

    this->prot->index_off();
    conf.insert(conf.begin()+insert_pos, _conf);
    conf[insert_pos].prot = this->prot;

    return insert_pos;
}

int Res::del_conf(const int& i_conf)
{
    if (i_conf<=0) return -1;   /* cannot delete backbone */
    if (i_conf>=conf.size()) return -1;

    this->prot->index_off();
    conf.erase(conf.begin()+i_conf);
    return 0;
}

string Res::resID()
{
    char _resID[10];
    sprintf(_resID, "%s_%c%04d%c", resName.c_str(), chainID, resSeq, iCode);
    int i_conf;
    for (i_conf=1;i_conf<conf.size();++i_conf) {
        if (fabs(conf[i_conf].H) > 1e-4) break;
        if (fabs(conf[i_conf].e) > 1e-4) break;
    }
    if (i_conf < conf.size()) {
        _resID[3] = conf[i_conf].confName[3];
    }

    string str(_resID);
    return str;
}

ostream& operator << (ostream& _file, Res& _res)
{
    _file << _res.resID();
    return _file;
}

inline vector <float> operator * (vector <float>& in_vec, float& in_num)
{
    vector <float> ret_vec(in_vec.size());
    for (int i=0;i<in_vec.size();++i) {
        ret_vec[i] = in_vec[i] * in_num;
    }
    return ret_vec;
}

void Res::make_sumcrg()
{
    int i_conf;
    print_sumcrg = false;
    for (i_conf=1; i_conf< conf.size(); ++i_conf) {
        if (fabs(conf[i_conf].H) > 1e-4) break;
        if (fabs(conf[i_conf].e) > 1e-4) break;
        if (fabs(conf[i_conf].netcrg) > 1e-4) break;
    }
    if (i_conf < conf.size()) print_sumcrg = true;

    if (conf.size() >= 2) {
        sum_H   = conf[1].occ_table * conf[1].H;
        sum_e   = conf[1].occ_table * conf[1].e;
        sum_crg = conf[1].occ_table * conf[1].netcrg;
    }
    for (int i_conf=2; i_conf< conf.size(); ++i_conf) {
        sum_H   = sum_2vec(conf[i_conf].occ_table * conf[i_conf].H, sum_H);
        sum_e   = sum_2vec(conf[i_conf].occ_table * conf[i_conf].e, sum_e);
        sum_crg = sum_2vec(conf[i_conf].occ_table * conf[i_conf].netcrg, sum_crg);
    }
}
