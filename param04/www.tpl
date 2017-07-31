# source for tpl?
CONFLIST WWW        WWWBK WWW+4 WWW+5 WWWDM

NATOM    WWWBK      0
NATOM    WWW+4      1
NATOM    WWW+5      1
NATOM	 WWW+6      1
NATOM    WWWDM      0

IATOM    WWW+4 W    0
IATOM    WWW+5 W    0
IATOM    WWW+6 W    0

ATOMNAME WWW+4   0  W
ATOMNAME WWW+5   0  W
ATOMNAME WWW+6   0  W

PROTON   WWW+4      0
pKA      WWW+4      0.0
ELECTRON WWW+4      0
EM       WWW_4      0
RXN      WWW+4      0

PROTON   WWW+5      0
PKA      WWW+5      0.0
ELECTRON WWW+5      1
EM       WWW+5      0
RXN      WWW+5      -186.75

PROTON   WWW+6      0
PKA      WWW+6      0.0
ELECTRON WWW+6      2
EM       WWW+6      0
RXN      WWW+6      0

# VDW_RAD and VDW_EPS also in 00always_needed.tpl
VDW_RAD  WWW+4 W     0
VDW_EPS  WWW+4 W     0
VDW_RAD  WWW+5 W     0
VDW_EPS  WWW+5 W     0

#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  WWW+4 W    ion
CONNECT  WWW+5 W    ion
CONNECT  WWW+6 W    ion

RADIUS   WWW    W   2.1  # source?

CHARGE   WWW+4 W    4.0
CHARGE   WWW+5 W    5.0
CHARGE   WWW+6 W    6.0

