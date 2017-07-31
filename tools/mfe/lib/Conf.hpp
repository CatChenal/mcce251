#include <string>
#include <vector>

using namespace std;

inline Conf::Conf()
{
}

inline Conf::Conf(const Res& _res)
{
    resName = _res.resName;
    chainID = _res.chainID;
    resSeq  = _res.resSeq;
    iCode   = _res.iCode;
}

inline Conf::Conf(const Conf& _conf)
{
    *this = _conf;
}

inline Conf::Conf(const Atom& _atom)
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

inline bool Conf::inRes(const Res& _res)
{
    if ( resName == _res.resName && chainID  == _res.chainID &&
        resSeq   == _res.resSeq && iCode   == _res.iCode )
        return true;
    else
        return false;
}

inline Conf& Conf::initialize_by_uniqID(const string& _uniqID)
{
    string msg = "wrong format for uniqID";
    if (_uniqID.size() != 14) throw msg;

    resName  = _uniqID.substr(0,3);
    confName = _uniqID.substr(0,5);
    chainID  = _uniqID[5];
    resSeq   = atoi(_uniqID.substr(6,4).c_str());
    iCode    = _uniqID[10];
    confID   = _uniqID.substr(11,3);
    uniqID   = _uniqID;
    
    return *this;
}

inline int Conf::ins_atom(const Atom& _atom)
{
    this->prot->index_off();
    atom.push_back(_atom);
    atom.back().prot = this->prot;
    
    return atom.size()-1;
}
