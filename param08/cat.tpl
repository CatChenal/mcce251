# Topology for catechol; 1,2-dihydroxybenzene; Chemical formula: C6H6O2
CONFLIST CAT        CATBK CAT01 CAT-1 CAT-2 CATDM

NATOM    CATBK      0
NATOM    CAT01      14
NATOM    CAT-1      13
NATOM    CAT-2      12
NATOM    CATDM      0

IATOM    CAT01  C1  0
IATOM    CAT01  O1  1
IATOM    CAT01  C2  2
IATOM    CAT01  O2  3
IATOM    CAT01  C3  4
IATOM    CAT01  C4  5
IATOM    CAT01  C5  6
IATOM    CAT01  C6  7
IATOM    CAT01 HO1  8
IATOM    CAT01  H2  9
IATOM    CAT01  H3  10
IATOM    CAT01  H4  11
IATOM    CAT01  H5  12
IATOM    CAT01 HO2  13

IATOM    CAT-1  C1  0
IATOM    CAT-1  O1  1
IATOM    CAT-1  C2  2
IATOM    CAT-1  O2  3
IATOM    CAT-1  C3  4
IATOM    CAT-1  C4  5
IATOM    CAT-1  C5  6
IATOM    CAT-1  C6  7
IATOM    CAT-1 HO1  8
IATOM    CAT-1  H2  9
IATOM    CAT-1  H3  10
IATOM    CAT-1  H4  11
IATOM    CAT-1  H5  12

IATOM    CAT-2  C1  0
IATOM    CAT-2  O1  1
IATOM    CAT-2  C2  2
IATOM    CAT-2  O2  3
IATOM    CAT-2  C3  4
IATOM    CAT-2  C4  5
IATOM    CAT-2  C5  6
IATOM    CAT-2  C6  7
IATOM    CAT-2  H2  8
IATOM    CAT-2  H3  9
IATOM    CAT-2  H4  10
IATOM    CAT-2  H5  11

ATOMNAME CAT01    0 C1
ATOMNAME CAT01    1 01
ATOMNAME CAT01    2 C2
ATOMNAME CAT01    3 O2
ATOMNAME CAT01    4 C3
ATOMNAME CAT01    5 C4
ATOMNAME CAT01    6 C5
ATOMNAME CAT01    7 C6
ATOMNAME CAT01    8 HO1
ATOMNAME CAT01    9 H2
ATOMNAME CAT01   10 H3
ATOMNAME CAT01   11 H4
ATOMNAME CAT01   12 H5
ATOMNAME CAT01   13 HO2

ATOMNAME CAT-1   0 C1
ATOMNAME CAT-1   1 O1
ATOMNAME CAT-1   2 C2
ATOMNAME CAT-1   3 O2
ATOMNAME CAT-1   4 C3
ATOMNAME CAT-1   5 C4
ATOMNAME CAT-1   6 C5
ATOMNAME CAT-1   7 C6
ATOMNAME CAT-1   8 HO1
ATOMNAME CAT-1   9 H2
ATOMNAME CAT-1  10 H3
ATOMNAME CAT-1  11 H4
ATOMNAME CAT-1  12 H5

ATOMNAME CAT-2   0 C1
ATOMNAME CAT-2   1 O1
ATOMNAME CAT-2   2 C2
ATOMNAME CAT-2   3 O2
ATOMNAME CAT-2   4 C3
ATOMNAME CAT-2   5 C4
ATOMNAME CAT-2   6 C5
ATOMNAME CAT-2   7 C6
ATOMNAME CAT-2   8 H2
ATOMNAME CAT-2   9 H3
ATOMNAME CAT-2  10 H4
ATOMNAME CAT-2  11 H5


#1.Basic Conformer Information: name, pka, em, rxn.

#23456789A123456789B123456789C
PROTON   CAT01      0
PROTON   CAT-1      -1
PROTON   CAT-2      -2
PROTON   CATDM      0

# Soln. pKa: pKa1 = 9.45; pKa2 = 12.8; CRC handbook of chemistry & physics
# https://pubchem.ncbi.nlm.nih.gov/compound/catechol#section=pKa

ELECTRON CAT01      0
ELECTRON CAT-1      0
ELECTRON CAT-2      0
EM       CAT01      0.0
EM       CAT-1      0.0
EM       CAT-2      0.0

PKA      CAT01      0
PKA      CAT-1      8.8
PKA      CAT-2      13.1
PKA      CATDM      0

RXN      CAT01      0
RXN      CAT-1      0
RXN      CAT-2      0


#2.Structure Connectivity

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  CAT01  C1  sp2       0      O1 0      C2 0      C6
CONNECT  CAT01  O1  sp3       0      C1 0     HO1
CONNECT  CAT01  C2  sp2       0      C1 0      O2 0      C3
CONNECT  CAT01  O2  sp3       0      C2 0     HO2
CONNECT  CAT01  C3  sp2       0      C2 0      C4 0      H2
CONNECT  CAT01  C4  sp2       0      C3 0      C5 0      H3
CONNECT  CAT01  C5  sp2       0      C4 0      C6 0      H4
CONNECT  CAT01  C6  sp2       0      C1 0      C5 0      H5
CONNECT  CAT01  HO1 s         0      O1
CONNECT  CAT01  H2  s         0      C3
CONNECT  CAT01  H3  s         0      C4
CONNECT  CAT01  H4  s         0      C5
CONNECT  CAT01  H5  s         0      C6
CONNECT  CAT01  HO2 s         0      O2

CONNECT  CAT-1  C1  sp2       0      O1 0      C2 0      C6
CONNECT  CAT-1  O1  sp3       0      C1 0     HO1
CONNECT  CAT-1  C2  sp2       0      C1 0      O2 0      C3
CONNECT  CAT-1  O2  sp2       0      C2
CONNECT  CAT-1  C3  sp2       0      C2 0      C4 0      H2
CONNECT  CAT-1  C4  sp2       0      C3 0      C5 0      H3
CONNECT  CAT-1  C5  sp2       0      C4 0      C6 0      H4
CONNECT  CAT-1  C6  sp2       0      C1 0      C5 0      H5
CONNECT  CAT-1  HO1 s         0      O1
CONNECT  CAT-1  H2  s         0      C3
CONNECT  CAT-1  H3  s         0      C4
CONNECT  CAT-1  H4  s         0      C5
CONNECT  CAT-1  H5  s         0      C6

CONNECT  CAT-2  C1  sp2       0      O1 0      C2 0      C6   
CONNECT  CAT-2  O1  sp2       0      C1    
CONNECT  CAT-2  C2  sp2       0      C1 0      O2 0      C3 
CONNECT  CAT-2  O2  sp2       0      C2 
CONNECT  CAT-2  C3  sp2       0      C2 0      C4 0      H2 
CONNECT  CAT-2  C4  sp2       0      C3 0      C5 0      H3 
CONNECT  CAT-2  C5  sp2       0      C4 0      C6 0      H4 
CONNECT  CAT-2  C6  sp2       0      C1 0      C5 0      H5 
CONNECT  CAT-2  H2  s         0      C3 
CONNECT  CAT-2  H3  s         0      C4 
CONNECT  CAT-2  H4  s         0      C5 
CONNECT  CAT-2  H5  s         0      C6 

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   CAT     C1 1.77
RADIUS   CAT     O1 1.45
RADIUS   CAT     C2 1.77
RADIUS   CAT     O2 1.45
RADIUS   CAT     C3 1.77
RADIUS   CAT     C4 1.77
RADIUS   CAT     C5 1.77
RADIUS   CAT     C6 1.77
RADIUS   CAT    HO1 0.60
RADIUS   CAT     H2 1.00
RADIUS   CAT     H3 1.00
RADIUS   CAT     H4 1.00
RADIUS   CAT     H5 1.00
RADIUS   CAT    HO2 0.60

CHARGE   CAT01   C1  0.269
CHARGE   CAT01   O1  -0.593
CHARGE   CAT01   C2  0.269
CHARGE   CAT01   O2  -0.593
CHARGE   CAT01   C3  -0.295
CHARGE   CAT01   C4  -0.084
CHARGE   CAT01   C5  -0.084
CHARGE   CAT01   C6  -0.295
CHARGE   CAT01  HO1  0.467
CHARGE   CAT01   H2  0.137
CHARGE   CAT01   H3  0.100
CHARGE   CAT01   H4  0.100
CHARGE   CAT01   H5  0.137
CHARGE   CAT01  HO2  0.467

CHARGE   CAT-1   C1  0.174
CHARGE   CAT-1   O1  -0.609
CHARGE   CAT-1   C2  0.362
CHARGE   CAT-1   O2  -0.788
CHARGE   CAT-1   C3  -0.412
CHARGE   CAT-1   C4  -0.053
CHARGE   CAT-1   C5  -0.195
CHARGE   CAT-1   C6  -0.285
CHARGE   CAT-1  HO1  0.444
CHARGE   CAT-1   H2  0.098
CHARGE   CAT-1   H3  0.069
CHARGE   CAT-1   H4  0.087
CHARGE   CAT-1   H5  0.108

CHARGE   CAT-2   C1  0.290
CHARGE   CAT-2   O1  -0.829
CHARGE   CAT-2   C2  0.290 
CHARGE   CAT-2   O2  -0.829
CHARGE   CAT-2   C3  -0.430
CHARGE   CAT-2   C4  -0.152
CHARGE   CAT-2   C5  -0.152
CHARGE   CAT-2   C6  -0.430
CHARGE   CAT-2   H2  0.067 
CHARGE   CAT-2   H3  0.054
CHARGE   CAT-2   H4  0.054
CHARGE   CAT-2   H5  0.067
