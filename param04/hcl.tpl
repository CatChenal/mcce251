# HCL renamed to HCL
#
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONFLIST HCL        HCLBK HCL-1 HCL01 HCLDM 

NATOM    HCLBK      0
NATOM    HCL-1      1
NATOM    HCL01      2
NATOM    HCLDM      0

IATOM    HCL-1 CL   0
ATOMNAME HCL-1    0 CL

IATOM    HCL01 CL   0
IATOM    HCL01  H   1
ATOMNAME HCL01    0 CL
ATOMNAME HCL01    1 H

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   HCL-1      0
PKA      HCL-1      0.0
ELECTRON HCL-1      0
EM       HCL-1      0.0
RXN      HCL-1      -20.374

PROTON   HCL01      1
PKA      HCL01      -6.0
ELECTRON HCL01      0
EM       HCL01      0.0
RXN      HCL01      0.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  HCL-1 CL   ion
CONNECT  HCL01 CL   ion
CONNECT  HCL01  H   s

#3.Atom Parameters: Partial Charges and Radii
RADIUS   HCL   CL   1.937 #Rashin, A. A., and Honig, B. (1985) Reevaluation of the Born Model of Ion Hydration, J Phys Chem-Us 89, 5588-5593.

CHARGE   HCL-1 CL   -1.0
CHARGE   HCL01 CL   0.0

RADIUS   HCL01  H   1.00
CHARGE   HCL01  H   0.0
