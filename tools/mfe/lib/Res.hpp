#include <string>
#include <vector>

using namespace std;

bool sameID(Conf &conf1, Conf &conf2);

inline Res::Res(const Atom& _atom)
{
    atomType    = _atom.atomType;
    chainID     = _atom.chainID;
    resName     = _atom.resName;
    resSeq      = _atom.resSeq;
    iCode       = _atom.iCode;
};

inline Res::Res(const Conf& _conf)
{
    chainID     = _conf.chainID;
    resName     = _conf.resName;
    resSeq      = _conf.resSeq;
    iCode       = _conf.iCode;
};

inline Conf& Res::asignbkb() {
    Conf bkb(*this);
    bkb.confName =  resName+"BK";
    bkb.confID   = "000";
    bkb.confHistory  = "BK________";
    bkb.altLoc   = ' ';
    
    if (conf.empty()) {
        ins_conf(bkb);
    }
    else if (!sameID(conf[0], bkb)) {
        conf.insert(conf.begin(), bkb);
    }
}

inline bool sameID(Conf &conf1, Conf &conf2) {
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

inline int Res::ins_conf(const Conf& _conf)
{
    this->prot->index_off();
    conf.push_back(_conf);
    conf.back().prot = this->prot;
    
    return conf.size()-1;
}

inline string Res::resID()
{
    char _resID[10];
    sprintf(_resID, "%s_%c%04d%c", resName.c_str(), chainID, resSeq, iCode);
    string str(_resID);
    return str;
}
