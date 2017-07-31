#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONFLIST _CL        _CLBK _CL-1 _CLDM 

NATOM    _CLBK      0
NATOM    _CL-1      1
NATOM    _CLDM      0

IATOM    _CL-1 CL   0

ATOMNAME _CL-1    0 CL  

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   _CLDM      0
PKA      _CLDM      0.0
ELECTRON _CLDM      0
EM       _CLDM      0.0
RXN      _CLDM      0

PROTON   _CL-1      0
PKA      _CL-1      0.0
ELECTRON _CL-1      0
EM       _CL-1      0.0

CHARGE   _CL-1 CL   -1.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  _CL-1 CL   ion


#3.Atom Parameters: Partial Charges and Radii

RXN      _CL-1      -20.36   #-20.323 : v251 public version
RADIUS   _CL   CL   1.937    # Rashin, A. A., and Honig, B. (1985) Reevaluation of the Born Model of Ion Hydration, J Phys Chem-Us 89, 5588-5593:
                             # ion	rad	dHhydra (kcal/mol)
                             # F-	1.423	-117.2
                             # Cl-	1.937	-86.1

# RXN      _CL-1      -21.925   # -21.925 goes with 1.80 rad
# RADIUS   _CL   CL   1.810 # Experimental ionic rad listed in Jensen & Jorgensen J. Chem. Theory. Comput. 2006 = Pauling radius

#-------|-----|----|--------------
#ParaNam|Res |Atom|Param/toggle
TRANS    _CL         t


