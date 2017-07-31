CONFLIST HOH        HOHBK HOH01 HOHDM #HOH-1 #HOH+1

NATOM    HOHBK      0
NATOM    HOH01      3
NATOM    HOHDM      0
NATOM    HOH-1      2
NATOM    HOH+1      4

IATOM    HOH01  O   0
IATOM    HOH01 1H   1
IATOM    HOH01 2H   2
IATOM    HOH-1  O   0
IATOM    HOH-1  H   1
IATOM    HOH+1  O   0
IATOM    HOH+1 1H   1
IATOM    HOH+1 2H   2
IATOM    HOH+1 3H   3

ATOMNAME HOH01    0  O
ATOMNAME HOH01    1 1H
ATOMNAME HOH01    2 2H
ATOMNAME HOH-1    0  O
ATOMNAME HOH-1    1  H
ATOMNAME HOH+1    0  O
ATOMNAME HOH+1    1 1H
ATOMNAME HOH+1    2 2H
ATOMNAME HOH+1    3 3H

#RXN      HOH01      0.0
RXN      HOH01      -0.796
RXN      HOH-1      0.0
RXN      HOH+1      0.0
RXN      HOHDM      0.0

PKA      HOH01      0.0
PKA      HOH-1      15.74
PKA      HOH+1      -1.74
PKA      HOHDM      0.0

PROTON   HOH01      0
PROTON   HOH-1      -1
PROTON   HOH+1      1
PROTON   HOHDM      0

ELECTRON HOH01      0
ELECTRON HOHDM      0

#23456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  HOH01  O   sp3       0    1H   0    2H
CONNECT  HOH01 1H   s         0     O
CONNECT  HOH01 2H   s         0     O

CONNECT  HOH-1  O   sp3       0     H
CONNECT  HOH-1  H   s         0     O

CONNECT  HOH+1  O   sp3       0    1H   0    2H   0    3H
CONNECT  HOH+1 1H   s         0     O
CONNECT  HOH+1 2H   s         0     O
CONNECT  HOH+1 3H   s         0     O

DONOR    HOH01 1H    O
DONOR    HOH01 2H    O
ACCEPTOR HOH01  O   1H

CHARGE   HOH01  O   -0.80
CHARGE   HOH01 1H   0.40
CHARGE   HOH01 2H   0.40

CHARGE   HOH-1  O   -1.20
CHARGE   HOH-1  H   0.20

CHARGE   HOH+1  O   -0.50
CHARGE   HOH+1 1H   0.50
CHARGE   HOH+1 2H   0.50
CHARGE   HOH+1 3H   0.50

RADIUS   HOH    O   1.60
RADIUS   HOH    H   1.00
RADIUS   HOH   1H   1.00
RADIUS   HOH   2H   1.00
RADIUS   HOH   3H   1.00
