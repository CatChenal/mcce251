#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "mcce.hpp"

using namespace std;

Conf::Conf()
{
}

Conf::Conf(const Res& _res)
{
    resName = _res.resName;
    chainID = _res.chainID;
    resSeq  = _res.resSeq;
    iCode   = _res.iCode;
}

Conf::Conf(const Conf& _conf)
{
    *this = _conf;
}

Conf::Conf(const Atom& _atom)
{
    resName = _atom.resName;
    chainID = _atom.chainID;
    resSeq  = _atom.resSeq;
    iCode   = _atom.iCode;

    confName = _atom.confName;
    confID  = _atom.confID;
    altLoc  = _atom.altLoc;
    confHistory = _atom.confHistory;
}

bool Conf::inRes(const Res& _res)
{
    if ( resName == _res.resName && chainID  == _res.chainID &&
        resSeq   == _res.resSeq && iCode   == _res.iCode )
        return true;
    else
        return false;
}

Conf& Conf::initialize_by_uniqID(const string& _uniqID)
{
    string msg = "Wrong format for uniqID (not 14 char long) ";
    if (_uniqID.size() != 14) {
        cerr << msg << _uniqID;
        throw msg;
    }
    resName  = _uniqID.substr(0,3);
    confName = _uniqID.substr(0,5);
    chainID  = _uniqID[5];
    resSeq   = atoi(_uniqID.substr(6,4).c_str());
    iCode    = _uniqID[10];
    confID   = _uniqID.substr(11,3);
    uniqID   = _uniqID;
    return *this;
}

int Conf::ins_atom(const Atom& _atom)
{
    this->prot->index_off();
    atom.push_back(_atom);
    atom.back().prot = this->prot;
    return atom.size()-1;
}

int Conf::find_katom(const string& _atomName)
{
    int i_atom;
    for (i_atom = 0; i_atom<atom.size(); i_atom++) {
        if (atom[i_atom].atomName == _atomName) break;
    }
    if (i_atom>= atom.size()) i_atom = -1;
    return i_atom;
}

string Conf::mk_uniqID()
{
    char temp_cstr[15];
    sprintf(temp_cstr, "%04d", resSeq);
    if (chainID == ' ') chainID = '_';
    if (iCode == ' ') iCode = '_';
    uniqID = confName + chainID + (string)temp_cstr + iCode + confID;
    return uniqID;
}
