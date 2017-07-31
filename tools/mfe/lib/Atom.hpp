#include <vector>
using namespace std;

inline Atom& Atom::read_pdbstr(const string& pdbstr)
{
    if (pdbstr.substr(0,6) != "ATOM  " && pdbstr.substr(0,6) != "HETATM") {
        switch_off();
        return *this;
    }
    
    switch_on();
    
    atomType    = (pdbstr.substr(0,6) == "ATOM  " ? ATOM : HETATM);     /* standard PDB format */
    atomSerial  = atoi(pdbstr.substr(6,5).c_str());                     /* standard PDB format */
    atomName    = pdbstr.substr(12,4);                                  /* standard PDB format */
    altLoc      = pdbstr[16];                                           /* standard PDB format */
    resName     = pdbstr.substr(17,3);                                  /* standard PDB format */
    chainID     = pdbstr[21];                                           /* standard PDB format */
    resSeq      = atoi(pdbstr.substr(22,4).c_str());                    /* standard PDB format */
    iCode       = pdbstr[26];                                           /* standard PDB format */
    for (int i=0;i<3;++i) r(i) = atof(pdbstr.substr(30+i*8,8).c_str());
    /*
    r(0)        = atof(pdbstr.substr(30,8).c_str());
    r(1)        = atof(pdbstr.substr(38,8).c_str());
    r(2)        = atof(pdbstr.substr(46,8).c_str());
    */
    
    rad         = atof(pdbstr.substr(55,7).c_str());
    crg         = atof(pdbstr.substr(68,6).c_str());

    confName    = pdbstr.substr(17,3) + pdbstr.substr(80,2);            /* MCCE */
    confID      = pdbstr.substr(27,3);                                  /* MCCE */
    confHistory = pdbstr.substr(80,10);                                 /* MCCE */
    
    return *this;
}

inline string& Atom::write_pdbstr(string& pdbstr)
{
    char cstring[160];
    pdbstr.erase();
    if (off()) return pdbstr;
    
    pdbstr = (atomType == ATOM ? "ATOM  " : "HETATM");
    sprintf(cstring, "%s%5d %4s%c%3s %c%04d%c%3s%8.3f%8.3f%8.3f %7.3f      %6.3f      %-10s",
    pdbstr.c_str(),
    atomSerial,
    atomName.c_str(),
    altLoc,
    resName.c_str(),
    chainID,
    resSeq,
    iCode,
    confID.c_str(),
    r[0],
    r[1],
    r[2],
    rad,
    crg,
    confHistory.c_str());
    
    pdbstr = cstring;
    return pdbstr;
}

inline bool Atom::inRes(const Res& _res)
{
    if ( resName == _res.resName && chainID  == _res.chainID &&
        resSeq   == _res.resSeq && iCode   == _res.iCode )
        return true;
    else
        return false;
}

inline bool Atom::inConf(const Conf& _conf)
{
    if ( resName == _conf.resName && chainID  == _conf.chainID &&
            resSeq   == _conf.resSeq && iCode   == _conf.iCode &&
            confName == _conf.confName && confID == _conf.confID)
        return true;
    else
        return false;
}

inline bool Atom::isbkb()
{
    if (confName.substr(3,2) == "BK")
        return true;
    else
        return false;
}

