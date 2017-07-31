#include <cstdlib>
#include <cstring>
#include <vector>
#include "mcce.hpp"

using namespace std;

Atom& Atom::read_pdbstr(const string& pdbstr)
{
    if (pdbstr.substr(0,6) != "ATOM  " && pdbstr.substr(0,6) != "HETATM") {
        switch_off();
        return *this;
    }

    switch_on();

    hetType    = (pdbstr.substr(0,6) == "ATOM  " ? ATOM : HETATM);     /* standard PDB format */
    atomSerial  = atoi(pdbstr.substr(6,5).c_str());                     /* standard PDB format */
    atomName    = pdbstr.substr(12,4);                                  /* standard PDB format */
    altLoc      = pdbstr[16];                                           /* standard PDB format */
    resName     = pdbstr.substr(17,3);                                  /* standard PDB format */
    chainID     = pdbstr[21];                                           /* standard PDB format */
    resSeq      = atoi(pdbstr.substr(22,4).c_str());                    /* standard PDB format */
    iCode       = pdbstr[26];                                           /* standard PDB format */
    for (int i=0;i<3;++i)
        r(i) = atof(pdbstr.substr(30+i*8,8).c_str());

    rad         = atof(pdbstr.substr(55,7).c_str());
    crg         = atof(pdbstr.substr(68,6).c_str());

    confName    = pdbstr.substr(17,3) + pdbstr.substr(80,2);            /* MCCE */
    confID      = pdbstr.substr(27,3);                                  /* MCCE */
    confHistory = pdbstr.substr(80,10);                                 /* MCCE */

    return *this;
}

Atom& Atom::read_orig_pdbstr(const string& pdbstr)
{
    if (pdbstr.substr(0,6) != "ATOM  " && pdbstr.substr(0,6) != "HETATM") {
        switch_off();
        return *this;
    }

    switch_on();

    hetType    = (pdbstr.substr(0,6) == "ATOM  " ? ATOM : HETATM);      /* standard PDB format */
    atomSerial  = atoi(pdbstr.substr(6,5).c_str());                     /* standard PDB format */
    atomName    = pdbstr.substr(12,4);                                  /* standard PDB format */
    altLoc      = pdbstr[16];                                           /* standard PDB format */
    resName     = pdbstr.substr(17,3);                                  /* standard PDB format */
    chainID     = pdbstr[21];                                           /* standard PDB format */
    resSeq      = atoi(pdbstr.substr(22,4).c_str());                    /* standard PDB format */
    iCode       = pdbstr[26];                                           /* standard PDB format */
    for (int i=0;i<3;++i)
        r(i) = atof(pdbstr.substr(30+i*8,8).c_str());                   /* standard PDB format */

    occ         = atof(pdbstr.substr(55,5).c_str());                    /* standard PDB format */
    bFactor     = atof(pdbstr.substr(61,5).c_str());                    /* standard PDB format */

    confName    = pdbstr.substr(17,3) + "  ";                           /* MCCE */
    confID      = pdbstr.substr(27,3);                                  /* MCCE */

    return *this;
}

string& Atom::write_pdbstr(string& pdbstr)
{
    char cstring[160];
    pdbstr.erase();
    if (off()) return pdbstr;

    pdbstr = (hetType == ATOM ? "ATOM  " : "HETATM");
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

Atom& Atom::read_demetristr(const string& pdbstr)
{
    if (pdbstr.substr(0,6) != "ATOM  " && pdbstr.substr(0,6) != "HETATM") {
        switch_off();
        return *this;
    }

    switch_on();
    //0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    //ATOM   AA    NT  NTR A 005A     11.078 -14.070 -15.860  -0.003
    hetType    = (pdbstr.substr(0,6) == "ATOM  " ? ATOM : HETATM);     /* standard PDB format */
    atomSerial  = atoi(pdbstr.substr(6,5).c_str());                     /* standard PDB format */
    atomName    = pdbstr.substr(12,4);                                  /* standard PDB format */
    altLoc      = pdbstr[16];                                           /* standard PDB format */
    resName     = pdbstr.substr(17,3);                                  /* standard PDB format */
    chainID     = pdbstr[21];                                           /* standard PDB format */
    resSeq      = atoi(pdbstr.substr(22,4).c_str());                    /* standard PDB format */
    iCode       = '_';                                                  /* standard PDB format */
    for (int i=0;i<3;++i)
        r(i) = atof(pdbstr.substr(30+i*8,8).c_str());

    crg         = atof(pdbstr.substr(54,8).c_str());                    /* demetri_out format */

    confName    = pdbstr.substr(17,3);                                  /* MCCE */
    char cstring[256];
    sprintf(cstring,"%03d",pdbstr[26]-'A'+1);
    confID      = cstring;                                              /* demetri */

    /* convert Res and Conf name */
    if (resName == "ASP") {resName = "ASP"; confName = "ASP01";}
    else if (resName == "AS-") {resName = "ASP"; confName = "ASP-1";}
    else if (resName == "CYS") {resName = "CYS"; confName = "CYS01";}
    else if (resName == "CY-") {resName = "CYS"; confName = "CYS-1";}
    else if (resName == "GLU") {resName = "GLU"; confName = "GLU01";}
    else if (resName == "GL-") {resName = "GLU"; confName = "GLU-1";}
    else if (resName == "TYR") {resName = "TYR"; confName = "TYR01";}
    else if (resName == "TY-") {resName = "TYR"; confName = "TYR-1";}
    else if (resName == "CTR") {resName = "CTR"; confName = "CTR01";}
    else if (resName == "CT-") {resName = "CTR"; confName = "CTR-1";}
    else if (resName == "ARG") {resName = "ARG"; confName = "ARG01";}
    else if (resName == "AR+") {resName = "ARG"; confName = "ARG+1";}
    else if (resName == "HIS") {resName = "HIS"; confName = "HIS01";}
    else if (resName == "HI+") {resName = "HIS"; confName = "HIS+1";}
    else if (resName == "LYS") {resName = "LYS"; confName = "LYS01";}
    else if (resName == "LY+") {resName = "LYS"; confName = "LYS+1";}
    else if (resName == "RSB") {resName = "RSB"; confName = "RSB01";}
    else if (resName == "RS+") {resName = "RSB"; confName = "RSB+1";}
    else if (resName == "NTR") {resName = "NTR"; confName = "NTR01";}
    else if (resName == "NT+") {resName = "NTR"; confName = "NTR+1";}
    else if (resName == "ALA") {resName = "ALA"; confName = "ALA01";}
    else if (resName == "ASN") {resName = "ASN"; confName = "ASN01";}
    else if (resName == "GLY") {resName = "GLY"; confName = "GLY01";}
    else if (resName == "GLN") {resName = "GLN"; confName = "GLN01";}
    else if (resName == "ILE") {resName = "ILE"; confName = "ILE01";}
    else if (resName == "LEU") {resName = "LEU"; confName = "LEU01";}
    else if (resName == "MET") {resName = "MET"; confName = "MET01";}
    else if (resName == "PHE") {resName = "PHE"; confName = "PHE01";}
    else if (resName == "PRO") {resName = "PRO"; confName = "PRO01";}
    else if (resName == "SER") {resName = "SER"; confName = "SER01";}
    else if (resName == "THR") {resName = "THR"; confName = "THR01";}
    else if (resName == "TRP") {resName = "TRP"; confName = "TRP01";}
    else if (resName == "VAL") {resName = "VAL"; confName = "VAL01";}

    /* RSB */
    if (atomName == " N0 ") atomName = " N  ";
    if (atomName == " H0 ") atomName = " HN ";
    if (atomName == " CA0") atomName = " CA ";
    if (atomName == " HA0") atomName = " HA ";
    if (atomName == " C0 ") atomName = " C  ";
    if (atomName == " O0 ") atomName = " O  ";

    if (confName == "RSB01") {
        if (atomName == " N  ") {switch_off();}
        if (atomName == " HN ") {switch_off();}
        if (atomName == " CA ") {switch_off();}
        if (atomName == " HA ") {switch_off();}
        if (atomName == " C  ") {switch_off();}
        if (atomName == " O  ") {switch_off();}
    }
    if (resName == "RSB") {
        if (atomName.substr(1,1) == "H") {switch_off();}
    }

    if (resName == "GLY" || atomName == " N  " || atomName == " HN " || atomName == " CA " || atomName == " HA " || atomName == " C  " || atomName == " O  ") {
        if (resName != "NTR" && resName != "CTR") {
            confID = "000";
            confName = resName + "BK";
        }
    }

    /* convert atom name */
    if (atomName == " HN ") atomName = " H  ";
    if (atomName == " HA1") atomName = "1HA ";
    if (atomName == " HA2") atomName = "2HA ";
    if (atomName == " HB1") atomName = "1HB ";
    if (atomName == " HB2") atomName = "2HB ";
    if (atomName == " HB3") atomName = "3HB ";
    if (atomName == " HG1") atomName = "1HG ";
    if (atomName == " HG2") atomName = "2HG ";
    if (resName != "PHE" && resName != "TYR" && resName != "TRP") {
        if (atomName == " HD1") atomName = "1HD ";
        if (atomName == " HD2") atomName = "2HD ";
    }
    if (resName != "PHE" && resName != "TYR" && resName != "TRP") {
        if (atomName == " HE1") atomName = "1HE ";
        if (atomName == " HE2") atomName = "2HE ";
        if (atomName == " HE3") atomName = "3HE ";
    }
    if (resName != "TRP") {
    if (atomName == " HZ1") atomName = "1HZ ";
    if (atomName == " HZ2") atomName = "2HZ ";
    if (atomName == " HZ3") atomName = "3HZ ";
    }
    if (atomName == " NT ") atomName = " N  ";
    if (atomName == " CAT") atomName = " CA ";
    if (atomName == " HT1") atomName = "1H  ";
    if (atomName == " HT2") atomName = "2H  ";
    if (atomName == " HT3") atomName = "3H  ";

    if (resName == "ASP") {
        if (atomName == " HT ") atomName = " HD1";
    }
    if (resName == "GLU") {
        if (atomName == " HT ") atomName = " HE1";
    }
    if (resName == "SER") {
        if (atomName == " HT ") atomName = " HG ";
    }
    if (resName == "THR") {
        if (atomName == " HT ") atomName = " HG1";
    }
    if (resName == "TYR") {
        if (atomName == " HT ") atomName = " HH ";
    }

    if (atomName == " HN ") atomName = " H  ";
    if (atomName == " HA1") atomName = "1HA ";
    if (atomName == " HA2") atomName = "2HA ";
    if (atomName == " HB1") atomName = "1HB ";
    if (atomName == " HB2") atomName = "2HB ";
    if (atomName == " HB3") atomName = "3HB ";
    if (atomName == " HG1") atomName = "1HG ";
    if (atomName == " HG2") atomName = "2HG ";

    if (resName == "MEM") {
        switch_off();
    }
    if (resName == "ARG") {
        if (atomName == " HH1" || atomName == " HH2" || atomName == " HE ") {
            switch_off();
        }
    }
    if (resName == "ASN") {
        if (atomName == "2HD ") {
            switch_off();
        }
    }
    if (resName == "GLN") {
        if (atomName == "2HE ") {
            switch_off();
        }
    }
    if (resName == "ILE") {
        if (atomName == "1HG " || atomName == "2HG " || atomName == "1HD ") {
            switch_off();
        }
    }
    if (resName == "LEU") {
        if (atomName == "1HD " || atomName == "2HD ") {
            switch_off();
        }
    }
    if (resName == "THR") {
        if (atomName == "1HG " || atomName == "2HG ") {
            switch_off();
        }
    }
    if (resName == "VAL") {
        if (atomName == "1HG " || atomName == "2HG ") {
            switch_off();
        }
    }
    if (confName == "LYS01") {
        if (atomName == "3HZ ") {
            switch_off();
        }
    }
    if (confName == "NTR01") {
        if (atomName == "3H  ") {
            switch_off();
        }
    }
    confHistory = confName.substr(3,2) + "________";

    return *this;
}

bool Atom::inRes(const Res& _res)
{
    if ( resName == _res.resName && chainID  == _res.chainID &&
        resSeq   == _res.resSeq && iCode   == _res.iCode )
    return true;
    else
        return false;
}

bool Atom::inConf(const Conf& _conf)
{
    if ( resName == _conf.resName && chainID  == _conf.chainID &&
        resSeq   == _conf.resSeq && iCode   == _conf.iCode &&
    confName == _conf.confName && confID == _conf.confID)
    return true;
    else
        return false;
}

bool Atom::isbkb()
{
    if (confName.substr(3,2) == "BK")
        return true;
    else
        return false;
}
