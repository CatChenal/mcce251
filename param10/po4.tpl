#-----------------------------------------------------------------------------
#     po4.tpl: Topology file for Phosphate:
#              phosphate, PO4(3-) (aka inorganic phosphate, or Pi)
#              hydrogen phosphate, HPO4(2-)
#              dihydrogen phosphate, H2PO4(1-)
#              phosphoric acid, H3PO4
#     Source: Giray Enkavi of Emad Tajkhorshid's lab (emad@life.uiuc.edu), 
#             except H3PO4, PO4(3-) (derived).
#     Created on 2012-04-17 (catalys)
#-----------------------------------------------------------------------------
CONFLIST PO4        PO4BK PO401 PO4-1 PO4-2 PO4-3 PO4DM

#23456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
#........|-----|----|---------|----|----|----|----|----|----|----|----|----|----|
NATOM    PO4BK      0
NATOM    PO401      8
NATOM    PO4-1      7
NATOM    PO4-2      6
NATOM    PO4-3      5
NATOM    PO4DM      0

IATOM    PO401  P   0
IATOM    PO401  O1  1
IATOM    PO401  HO1 2
IATOM    PO401  O2  3
IATOM    PO401  HO2 4
IATOM    PO401  O3  5
IATOM    PO401  HO3 6
IATOM    PO401  O4  7

IATOM    PO4-1  P   0
IATOM    PO4-1  O1  1
IATOM    PO4-1  HO1 2
IATOM    PO4-1  O2  3
IATOM    PO4-1  HO2 4
IATOM    PO4-1  O3  5
IATOM    PO4-1  O4  6

IATOM    PO4-2  P   0
IATOM    PO4-2  O1  1
IATOM    PO4-2  HO1 2
IATOM    PO4-2  O2  3
IATOM    PO4-2  O3  4
IATOM    PO4-2  O4  5

IATOM    PO4-3  P   0
IATOM    PO4-3  O1  1
IATOM    PO4-3  O2  2
IATOM    PO4-3  O3  3
IATOM    PO4-3  O4  4

ATOMNAME PO401    0  P
ATOMNAME PO401    1  O1 
ATOMNAME PO401    2  HO1 
ATOMNAME PO401    3  O2 
ATOMNAME PO401    4  HO2 
ATOMNAME PO401    5  O3 
ATOMNAME PO401    6  HO3
ATOMNAME PO401    7  O4 

ATOMNAME PO4-1    0  P
ATOMNAME PO4-1    1  O1 
ATOMNAME PO4-1    2  HO1 
ATOMNAME PO4-1    3  O2 
ATOMNAME PO4-1    4  HO2 
ATOMNAME PO4-1    5  O3 
ATOMNAME PO4-1    6  O4 

ATOMNAME PO4-2    0  P
ATOMNAME PO4-2    1  O1 
ATOMNAME PO4-2    2  HO1 
ATOMNAME PO4-2    3  O2 
ATOMNAME PO4-2    4  O3 
ATOMNAME PO4-2    5  O4 

ATOMNAME PO4-3    0  P
ATOMNAME PO4-3    1  O1 
ATOMNAME PO4-3    2  O2 
ATOMNAME PO4-3    3  O3 
ATOMNAME PO4-3    4  O4 

# 1. Basic Conformer Information:
# Number of protons and electrons, pKa, Em and Reaction Field Energy (RXN)

#23456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
# 'PROTON' means charge here:
PROTON   PO401      0
PROTON   PO4-1      -1
PROTON   PO4-2      -2
PROTON   PO4-3      -3

# Solution pKa1=2.12, pKa2=7.21, and pKa3=12.67 : CRC Handbook of Chemistry and Physics
PKA      PO401      0.0
PKA      PO4-1      2.12
PKA      PO4-2      7.21
PKA      PO4-3      12.67

ELECTRON PO401      0
ELECTRON PO4-1      0
ELECTRON PO4-2      0
ELECTRON PO4-3      0

EM       PO401      0.0
EM       PO4-1      0.0
EM       PO4-2      0.0
EM       PO4-3      0.0

RXN      PO401      -2.600
RXN      PO4-1      -8.366
RXN      PO4-2      -26.362
RXN      PO4-3      -56.456

# H3PO4
#23456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  PO401  P   sp3       0     O1  0     O2  0     O3  0     O4
CONNECT  PO401  O1  sp3       0     P   0     HO1
CONNECT  PO401  HO1 s         0     O1
CONNECT  PO401  O2  sp3       0     P   0     HO2
CONNECT  PO401  HO2 s         0     O2
CONNECT  PO401  O3  sp3       0     P   0     HO3
CONNECT  PO401  HO3 s         0     O3
CONNECT  PO401  O4  sp2       0     P

# H2PO4-
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  PO4-1  P   sp3       0     O1  0     O2  0     O3  0     O4
CONNECT  PO4-1  O1  sp3       0     P   0     HO1 
CONNECT  PO4-1  HO1 s         0     O1 
CONNECT  PO4-1  O2  sp3       0     P   0     HO2
CONNECT  PO4-1  HO2 s         0     O2
CONNECT  PO4-1  O3  sp3       0     P
CONNECT  PO4-1  O4  sp2       0     P

# HPO42-
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  PO4-2  P   sp3       0     O1  0     O2  0     O3  0     O4
CONNECT  PO4-2  O1  sp3       0     P   0     HO1
CONNECT  PO4-2  HO1 s         0     O1
CONNECT  PO4-2  O2  sp3       0     P
CONNECT  PO4-2  O3  sp3       0     P
CONNECT  PO4-2  O4  s         0     P

# PO43- : phosphate
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  PO4-3  P   sp3       0     O1  0     O2  0     O3  0     O4
CONNECT  PO4-3  O1  sp3       0     P
CONNECT  PO4-3  O2  sp3       0     P
CONNECT  PO4-3  O3  sp3       0     P
CONNECT  PO4-3  O4  sp2       0     P

# Acceptor/donor section:
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#-------|-----|----|----|----|----|----|---------|---------|---------|----
#        CONF  ATOM ATOM ATOM ATOM      phi0(min)  n_fold   Amplitude(barrier,kcal/mol)     
DONOR    PO401 HO1   O1 
DONOR    PO401 HO2   O2 
DONOR    PO401 HO3   O3 
ACCEPTOR PO4-1 O3 
DONOR    PO4-1 HO1   O1 
DONOR    PO4-1 HO2   O2 
ACCEPTOR PO4-2 O2 
ACCEPTOR PO4-2 O3 
DONOR    PO4-2 HO1   O1 
ACCEPTOR PO4-3 O1 
ACCEPTOR PO4-3 O2 
ACCEPTOR PO4-3 O3 

# 3. Atom Parameters: 
#    Partial Charges: Giray Enkavi (CHARMM)
#    Radii          : "Bondi, J.Phys.Chem., 68, 441, 1964."
CHARGE   PO401  P     1.92
CHARGE   PO401  O1   -0.68
CHARGE   PO401  HO1   0.34
CHARGE   PO401  O2   -0.68
CHARGE   PO401  HO2   0.34
CHARGE   PO401  O3   -0.68
CHARGE   PO401  HO3   0.34
CHARGE   PO401  O4   -0.90

CHARGE   PO4-1  P     1.48
CHARGE   PO4-1  O1   -0.68
CHARGE   PO4-1  HO1   0.34
CHARGE   PO4-1  O2   -0.68
CHARGE   PO4-1  HO2   0.34
CHARGE   PO4-1  O3   -0.90
CHARGE   PO4-1  O4   -0.90

CHARGE   PO4-2  P     1.04
CHARGE   PO4-2  O1   -0.68
CHARGE   PO4-2  HO1   0.34
CHARGE   PO4-2  O2   -0.90
CHARGE   PO4-2  O3   -0.90
CHARGE   PO4-2  O4   -0.90

CHARGE   PO4-3  P     0.60
CHARGE   PO4-3  O1   -0.90
CHARGE   PO4-3  O2   -0.90
CHARGE   PO4-3  O3   -0.90
CHARGE   PO4-3  O4   -0.90

RADIUS   PO401  P     1.90
RADIUS   PO401  O1    1.40
RADIUS   PO401  O2    1.40
RADIUS   PO401  O3    1.40
RADIUS   PO401  O4    1.40
RADIUS   PO401  HO1   1.00
RADIUS   PO401  HO2   1.00
RADIUS   PO401  HO3   1.0

RADIUS   PO4-1  P     1.90
RADIUS   PO4-1  O1    1.40
RADIUS   PO4-1  O2    1.40
RADIUS   PO4-1  O3    1.40
RADIUS   PO4-1  O4    1.40
RADIUS   PO4-1  HO1   1.00
RADIUS   PO4-1  HO2   1.00

RADIUS   PO4-2  P     1.90
RADIUS   PO4-2  O1    1.40
RADIUS   PO4-2  O2    1.40
RADIUS   PO4-2  O3    1.40
RADIUS   PO4-2  O4    1.40
RADIUS   PO4-2  HO1   1.00

RADIUS   PO4-3  P     1.90
RADIUS   PO4-3  O1    1.40
RADIUS   PO4-3  O2    1.40
RADIUS   PO4-3  O3    1.40
RADIUS   PO4-3  O4    1.40

# VDW params: from debug.log:
VDW_RAD  PO401  P     1.91
VDW_EPS  PO401  P     0.11
VDW_RAD  PO401  O1    1.66
VDW_EPS  PO401  O1    0.21
VDW_RAD  PO401  O2    1.66
VDW_EPS  PO401  O2    0.21
VDW_RAD  PO401  O3    1.66
VDW_EPS  PO401  O3    0.21
VDW_RAD  PO401  O4    1.66
VDW_EPS  PO401  O4    0.21
VDW_RAD  PO401  HO1   1.49
VDW_EPS  PO401  HO1   0.02
VDW_RAD  PO401  HO2   1.49
VDW_EPS  PO401  HO2   0.02
VDW_RAD  PO401  HO3   1.49
VDW_EPS  PO401  HO3   0.02

VDW_RAD  PO4-1  P     1.91
VDW_EPS  PO4-1  P     0.11
VDW_RAD  PO4-1  O1    1.66
VDW_EPS  PO4-1  O1    0.21
VDW_RAD  PO4-1  HO1   1.49
VDW_EPS  PO4-1  HO1   0.02
VDW_RAD  PO4-1  O2    1.66
VDW_EPS  PO4-1  O2    0.21
VDW_RAD  PO4-1  HO2   1.49
VDW_EPS  PO4-1  HO2   0.02
VDW_RAD  PO4-1  O3    1.66
VDW_EPS  PO4-1  O3    0.21
VDW_RAD  PO4-1  O4    1.66
VDW_EPS  PO4-1  O4    0.21

VDW_RAD  PO4-2  P     1.91
VDW_EPS  PO4-2  P     0.11
VDW_RAD  PO4-2  O1    1.66
VDW_EPS  PO4-2  O1    0.21
VDW_RAD  PO4-2  HO1   1.49
VDW_EPS  PO4-2  HO1   0.02
VDW_RAD  PO4-2  O2    1.66
VDW_EPS  PO4-2  O2    0.21
VDW_RAD  PO4-2  O3    1.66
VDW_EPS  PO4-2  O3    0.21
VDW_RAD  PO4-2  O4    1.66
VDW_EPS  PO4-2  O4    0.21

VDW_RAD  PO4-3  P     1.91
VDW_EPS  PO4-3  P     0.11
VDW_RAD  PO4-3  O1    1.66
VDW_EPS  PO4-3  O1    0.21
VDW_RAD  PO4-3  O2    1.66
VDW_EPS  PO4-3  O2    0.21
VDW_RAD  PO4-3  O3    1.66
VDW_EPS  PO4-3  O3    0.21
VDW_RAD  PO4-3  O4    1.66
VDW_EPS  PO4-3  O4    0.21

TORSION  PO401  HO1  O1   P    O2   f        0.0         3      0.00
TORSION  PO401  HO2  O2   P    O1   f        0.0         3      0.00
TORSION  PO401  HO3  O3   P    O1   f        0.0         3      0.00
TORSION  PO4-2  HO1  O1   P    O2   f        0.0         3      0.00
TORSION  PO4-1  HO1  O1   P    O2   f        0.0         3      0.00
TORSION  PO4-1  HO2  O2   P    O1   f        0.0         3      0.00

# SWAP RULE
#ROT_SWAP HIS   0     ND1- CD2  CE1- NE2
#23456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
#--------|-----|----|----------|---------|---------|---------|---------|---------|
ROT_SWAP PO4   0     O4 - O1
ROT_SWAP PO4   1     O4 - O2
ROT_SWAP PO4   2     O4 - O3
ROT_SWAP PO4   3     O3 - O2
ROT_SWAP PO4   4     O2 - O1
