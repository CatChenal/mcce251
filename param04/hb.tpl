#23456789A123456789B123456789C123456789D123456789
HDONOR   ASP01       HD1
HDONOR   ASP02       HD2

HDONOR   ARG01      1HH1 2HH1 1HH2 2HH2
HDONOR   ARG02      1HH1 1HH2 2HH2  HE
HDONOR   ARG03      1HH1 2HH1 1HH2  HE
HDONOR   ARG+1      1HH1 2HH1 1HH2 2HH2  HE

HDONOR   ASN01      1HD2 2HD2

HDONOR   CYS01       HG
 
HDONOR   GLN01      1HE2 2HE2

HDONOR   GLU01       HE1
HDONOR   GLU02       HE2

HDONOR   HIS01       HE2
HDONOR   HIS02       HD1
HDONOR   HIS+1       HD1  HE2

HDONOR   SER01       HG 

HDONOR   TRP01       HE1

HDONOR   THR01       HG1

HDONOR   TYR01       HH

HDONOR   LYS01      1HZ  2HZ
HDONOR   LYS+1      1HZ  2HZ  3HZ

HDONOR   PAA01       H1A
HDONOR   PAA02       H2A

HDONOR   PDD01       H1D
HDONOR   PDD02       H2A
################################################
# hydrogen bond acceptor #######################
################################################

HACCEPT  ASN01       OD1 

HACCEPT  ASP01       OD2
HACCEPT  ASP02       OD1
HACCEPT  ASP-1       OD1  OD2

HACCEPT  GLN01       OE1

HACCEPT  GLU01       OE2
HACCEPT  GLU02       OE1
HACCEPT  GLU-1       OE1  OE2

HACCEPT  TYR-1       OH

HACCEPT  PAA01       O2A
HACCEPT  PAA02       O1A
HACCEPT  PAA-1       O1A  O2A

HACCEPT  PDD01       O2D
HACCEPT  PDD02       O1D
HACCEPT  PDD-1       O1D  O2D

HACCEPT  MET01       SD


# not often used amino acid 
HDONOR   UbQP1       HO2
HDONOR   UbQP2       HO5
HDONOR   UbQS1       HO2
HDONOR   UbQS2       HO5
HDONOR   HOH01      1H   2H          	# water molecule
HDONOR   HOH-1      1H               	# water molecule
HDONOR   HOH+1      1H   2H   3H        # water molecule

HDONOR   PTR01      HO2P HO3P
HDONOR   PTR-1      HO2P
HDONOR   TPO01      HO2P HO3P
HDONOR   TPO-1      HO2P
HDONOR   TYL01       HH

HACCEPT  HOH01       O			# water molecule
HACCEPT  HOH-1       O			# water molecule
HACCEPT  HOH+1       O			# water molecule
HACCEPT  PTR-1       O2P
HACCEPT  PTR-2       O2P  O3P
HACCEPT  TPO-1       O2P
HACCEPT  TPO-2       O2P  O3P
HACCEPT  TYL-1       OH

# for rsb
HDONOR   RSB+1       HZ
HACCEPT  RSB01       NZ
HACCEPT  RSB+1       NZ


# for HE3
#HDONOR   HE3+O       HD1 
#HACCEPT  HE3+O	     ND1
#HDONOR   HE30W       HD1 1H   2H
#HACCEPT  HE30W       ND1  O
#HDONOR   HE3+W       HD1 1H   2H
#HACCEPT  HE3+W       ND1  O
#HDONOR   HE30H       HD1  H
#HACCEPT  HE30W       ND1  O
#HDONOR   HE3+H       HD1  H
#HACCEPT  HE3+H       ND1  O 


# for HIS 
HACCEPT   HIS01      ND1  NE2
HACCEPT   HIS02      ND1  NE2
HACCEPT   HIS+1      ND1  NE2
