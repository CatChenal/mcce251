CONFLIST WII        WIIBK WII+4 WII+5 WIIDM

NATOM    WIIBK      0
NATOM    WII+4      47
NATOM    WII+5      47
NATOM    WIIDM      0

IATOM    WII+4 BR1  0
IATOM    WII+4 C01  1
IATOM    WII+4 C02  2
IATOM    WII+4 C03  3
IATOM    WII+4 C04  4
IATOM    WII+4 C05  5
IATOM    WII+4 C06  6
IATOM    WII+4 C07  7
IATOM    WII+4 H01  8
IATOM    WII+4 H02  9
IATOM    WII+4 H03  10
IATOM    WII+4 H04  11
IATOM    WII+4 H05  12
IATOM    WII+4 C08  13
IATOM    WII+4 O01  14
IATOM    WII+4 C09  15
IATOM    WII+4 C10  16
IATOM    WII+4 C11  17
IATOM    WII+4 O02  18
IATOM    WII+4 N01  19
IATOM    WII+4 N02  20
IATOM    WII+4 C12  21
IATOM    WII+4 C13  22
IATOM    WII+4 C14  23
IATOM    WII+4 H06  24
IATOM    WII+4 H07  25
IATOM    WII+4 H08  26
IATOM    WII+4 C15  27
IATOM    WII+4 H09  28
IATOM    WII+4 H10  29
IATOM    WII+4 H11  30
IATOM    WII+4 C16  31
IATOM    WII+4 H12  32
IATOM    WII+4 H13  33
IATOM    WII+4 H14  34
IATOM    WII+4 C17  35
IATOM    WII+4 H15  36
IATOM    WII+4 H16  37
IATOM    WII+4 H17  38
IATOM    WII+4 C18  39
IATOM    WII+4 H18  40
IATOM    WII+4 H19  41
IATOM    WII+4 H20  42
IATOM    WII+4 C19  43
IATOM    WII+4 H21  44
IATOM    WII+4 H22  45
IATOM    WII+4 H23  46
ATOMNAME WII+4   0  BR1
ATOMNAME WII+4   1  C01
ATOMNAME WII+4   2  C02
ATOMNAME WII+4   3  C03
ATOMNAME WII+4   4  C04
ATOMNAME WII+4   5  C05
ATOMNAME WII+4   6  C06
ATOMNAME WII+4   7  C07
ATOMNAME WII+4   8  H01
ATOMNAME WII+4   9  H02
ATOMNAME WII+4   10  H03
ATOMNAME WII+4   11  H04
ATOMNAME WII+4   12  H05
ATOMNAME WII+4   13  C08
ATOMNAME WII+4   14  O01
ATOMNAME WII+4   15  C09
ATOMNAME WII+4   16  C10
ATOMNAME WII+4   17  C11
ATOMNAME WII+4   18  O02
ATOMNAME WII+4   19  N01
ATOMNAME WII+4   20  N02
ATOMNAME WII+4   21  C12
ATOMNAME WII+4   22  C13
ATOMNAME WII+4   23  C14
ATOMNAME WII+4   24  H06
ATOMNAME WII+4   25  H07
ATOMNAME WII+4   26  H08
ATOMNAME WII+4   27  C15
ATOMNAME WII+4   28  H09
ATOMNAME WII+4   29  H10
ATOMNAME WII+4   30  H11
ATOMNAME WII+4   31  C16
ATOMNAME WII+4   32  H12
ATOMNAME WII+4   33  H13
ATOMNAME WII+4   34  H14
ATOMNAME WII+4   35  C17
ATOMNAME WII+4   36  H15
ATOMNAME WII+4   37  H16
ATOMNAME WII+4   38  H17
ATOMNAME WII+4   39  C18
ATOMNAME WII+4   40  H18
ATOMNAME WII+4   41  H19
ATOMNAME WII+4   42  H20
ATOMNAME WII+4   43  C19
ATOMNAME WII+4   44  H21
ATOMNAME WII+4   45  H22
ATOMNAME WII+4   46  H23
CONNECT  WII+4  BR1 ion   
CONNECT  WII+4  C01 ion   0   C02   
CONNECT  WII+4  C02 ion   0   C01   0   C03   0   C07   0   H01   0   H05   
CONNECT  WII+4  C03 ion   0   C02   0   C04   0   H01   0   H02   
CONNECT  WII+4  C04 ion   0   C03   0   C05   0   H01   0   H02   0   H03   
CONNECT  WII+4  C05 ion   0   C04   0   C06   0   H02   0   H03   0   H04   
CONNECT  WII+4  C06 ion   0   C05   0   C07   0   H03   0   H04   0   H05   
CONNECT  WII+4  C07 ion   0   C02   0   C06   0   H04   0   H05   
CONNECT  WII+4  H01 ion   0   C02   0   C03   0   C04   
CONNECT  WII+4  H02 ion   0   C03   0   C04   0   C05   
CONNECT  WII+4  H03 ion   0   C04   0   C05   0   C06   
CONNECT  WII+4  H04 ion   0   C05   0   C06   0   C07   
CONNECT  WII+4  H05 ion   0   C02   0   C06   0   C07   
CONNECT  WII+4  C08 ion   0   O01   
CONNECT  WII+4  O01 ion   0   C08   
CONNECT  WII+4  C09 ion   0   N01   
CONNECT  WII+4  C10 ion   0   N02   
CONNECT  WII+4  C11 ion   0   O02   
CONNECT  WII+4  O02 ion   0   C11   
CONNECT  WII+4  N01 ion   0   C09   0   C13   
CONNECT  WII+4  N02 ion   0   C10   0   C12   
CONNECT  WII+4  C12 ion   0   N02   0   C17   0   H15   0   H16   0   H17   0   C18   0   H18   0   H19   0   H20   0   C19   0   H21   0   H22   0   H23   
CONNECT  WII+4  C13 ion   0   N01   0   C14   0   H06   0   H07   0   H08   0   C15   0   H09   0   H10   0   H11   0   C16   0   H12   0   H13   0   H14   
CONNECT  WII+4  C14 ion   0   C13   0   H06   0   H07   0   H08   
CONNECT  WII+4  H06 ion   0   C13   0   C14   0   H07   0   H08   
CONNECT  WII+4  H07 ion   0   C13   0   C14   0   H06   0   H08   
CONNECT  WII+4  H08 ion   0   C13   0   C14   0   H06   0   H07   
CONNECT  WII+4  C15 ion   0   C13   0   H09   0   H10   0   H11   
CONNECT  WII+4  H09 ion   0   C13   0   C15   0   H10   0   H11   
CONNECT  WII+4  H10 ion   0   C13   0   C15   0   H09   0   H11   
CONNECT  WII+4  H11 ion   0   C13   0   C15   0   H09   0   H10   
CONNECT  WII+4  C16 ion   0   C13   0   H12   0   H13   0   H14   
CONNECT  WII+4  H12 ion   0   C13   0   C16   0   H13   0   H14   
CONNECT  WII+4  H13 ion   0   C13   0   C16   0   H12   0   H14   
CONNECT  WII+4  H14 ion   0   C13   0   C16   0   H12   0   H13   
CONNECT  WII+4  C17 ion   0   C12   0   H15   0   H16   0   H17   
CONNECT  WII+4  H15 ion   0   C12   0   C17   0   H16   0   H17   
CONNECT  WII+4  H16 ion   0   C12   0   C17   0   H15   0   H17   
CONNECT  WII+4  H17 ion   0   C12   0   C17   0   H15   0   H16   
CONNECT  WII+4  C18 ion   0   C12   0   H18   0   H19   0   H20   
CONNECT  WII+4  H18 ion   0   C12   0   C18   0   H19   0   H20   
CONNECT  WII+4  H19 ion   0   C12   0   C18   0   H18   0   H20   
CONNECT  WII+4  H20 ion   0   C12   0   C18   0   H18   0   H19   
CONNECT  WII+4  C19 ion   0   C12   0   H21   0   H22   0   H23   
CONNECT  WII+4  H21 ion   0   C12   0   C19   0   H22   0   H23   
CONNECT  WII+4  H22 ion   0   C12   0   C19   0   H21   0   H23   
CONNECT  WII+4  H23 ion   0   C12   0   C19   0   H21   0   H22   
CHARGE   WII+4 BR1  -0.815954
CHARGE   WII+4 C01  -1.458085
CHARGE   WII+4 C02  0.647091
CHARGE   WII+4 C03  -0.291433
CHARGE   WII+4 C04  -0.157535
CHARGE   WII+4 C05  -0.128369
CHARGE   WII+4 C06  -0.211740
CHARGE   WII+4 C07  -0.230923
CHARGE   WII+4 H01  0.180115
CHARGE   WII+4 H02  0.166704
CHARGE   WII+4 H03  0.156408
CHARGE   WII+4 H04  0.174035
CHARGE   WII+4 H05  0.190743
CHARGE   WII+4 C08  -0.607095
CHARGE   WII+4 O01  -0.046423
CHARGE   WII+4 C09  -1.165039
CHARGE   WII+4 C10  -1.176434
CHARGE   WII+4 C11  -0.599222
CHARGE   WII+4 O02  -0.049721
CHARGE   WII+4 N01  0.456308
CHARGE   WII+4 N02  0.465519
CHARGE   WII+4 C12  0.618436
CHARGE   WII+4 C13  0.610123
CHARGE   WII+4 C14  -0.489082
CHARGE   WII+4 H06  0.140677
CHARGE   WII+4 H07  0.118418
CHARGE   WII+4 H08  0.132041
CHARGE   WII+4 C15  -0.568143
CHARGE   WII+4 H09  0.139505
CHARGE   WII+4 H10  0.168689
CHARGE   WII+4 H11  0.132766
CHARGE   WII+4 C16  -0.634117
CHARGE   WII+4 H12  0.166656
CHARGE   WII+4 H13  0.182443
CHARGE   WII+4 H14  0.152328
CHARGE   WII+4 C17  -0.585959
CHARGE   WII+4 H15  0.171769
CHARGE   WII+4 H16  0.140562
CHARGE   WII+4 H17  0.143004
CHARGE   WII+4 C18  -0.586004
CHARGE   WII+4 H18  0.147080
CHARGE   WII+4 H19  0.165610
CHARGE   WII+4 H20  0.155126
CHARGE   WII+4 C19  -0.569863
CHARGE   WII+4 H21  0.168688
CHARGE   WII+4 H22  0.143545
CHARGE   WII+4 H23  0.136755
RADIUS   WII+4  BR1   1.85
RADIUS   WII+4  C01   1.70
RADIUS   WII+4  C02   1.70
RADIUS   WII+4  C03   1.70
RADIUS   WII+4  C04   1.70
RADIUS   WII+4  C05   1.70
RADIUS   WII+4  C06   1.70
RADIUS   WII+4  C07   1.70
RADIUS   WII+4  H01   1.20
RADIUS   WII+4  H02   1.20
RADIUS   WII+4  H03   1.20
RADIUS   WII+4  H04   1.20
RADIUS   WII+4  H05   1.20
RADIUS   WII+4  C08   1.70
RADIUS   WII+4  O01   1.52
RADIUS   WII+4  C09   1.70
RADIUS   WII+4  C10   1.70
RADIUS   WII+4  C11   1.70
RADIUS   WII+4  O02   1.52
RADIUS   WII+4  N01   1.55
RADIUS   WII+4  N02   1.55
RADIUS   WII+4  C12   1.70
RADIUS   WII+4  C13   1.70
RADIUS   WII+4  C14   1.70
RADIUS   WII+4  H06   1.20
RADIUS   WII+4  H07   1.20
RADIUS   WII+4  H08   1.20
RADIUS   WII+4  C15   1.70
RADIUS   WII+4  H09   1.20
RADIUS   WII+4  H10   1.20
RADIUS   WII+4  H11   1.20
RADIUS   WII+4  C16   1.70
RADIUS   WII+4  H12   1.20
RADIUS   WII+4  H13   1.20
RADIUS   WII+4  H14   1.20
RADIUS   WII+4  C17   1.70
RADIUS   WII+4  H15   1.20
RADIUS   WII+4  H16   1.20
RADIUS   WII+4  H17   1.20
RADIUS   WII+4  C18   1.70
RADIUS   WII+4  H18   1.20
RADIUS   WII+4  H19   1.20
RADIUS   WII+4  H20   1.20
RADIUS   WII+4  C19   1.70
RADIUS   WII+4  H21   1.20
RADIUS   WII+4  H22   1.20
RADIUS   WII+4  H23   1.20
IATOM    WII+5 BR1  0
IATOM    WII+5 C01  1
IATOM    WII+5 C02  2
IATOM    WII+5 C03  3
IATOM    WII+5 C04  4
IATOM    WII+5 C05  5
IATOM    WII+5 C06  6
IATOM    WII+5 C07  7
IATOM    WII+5 H01  8
IATOM    WII+5 H02  9
IATOM    WII+5 H03  10
IATOM    WII+5 H04  11
IATOM    WII+5 H05  12
IATOM    WII+5 C08  13
IATOM    WII+5 O01  14
IATOM    WII+5 C09  15
IATOM    WII+5 C10  16
IATOM    WII+5 C11  17
IATOM    WII+5 O02  18
IATOM    WII+5 N01  19
IATOM    WII+5 N02  20
IATOM    WII+5 C12  21
IATOM    WII+5 C13  22
IATOM    WII+5 C14  23
IATOM    WII+5 H06  24
IATOM    WII+5 H07  25
IATOM    WII+5 H08  26
IATOM    WII+5 C15  27
IATOM    WII+5 H09  28
IATOM    WII+5 H10  29
IATOM    WII+5 H11  30
IATOM    WII+5 C16  31
IATOM    WII+5 H12  32
IATOM    WII+5 H13  33
IATOM    WII+5 H14  34
IATOM    WII+5 C17  35
IATOM    WII+5 H15  36
IATOM    WII+5 H16  37
IATOM    WII+5 H17  38
IATOM    WII+5 C18  39
IATOM    WII+5 H18  40
IATOM    WII+5 H19  41
IATOM    WII+5 H20  42
IATOM    WII+5 C19  43
IATOM    WII+5 H21  44
IATOM    WII+5 H22  45
IATOM    WII+5 H23  46

ATOMNAME WII+5   0  BR1
ATOMNAME WII+5   1  C01
ATOMNAME WII+5   2  C02
ATOMNAME WII+5   3  C03
ATOMNAME WII+5   4  C04
ATOMNAME WII+5   5  C05
ATOMNAME WII+5   6  C06
ATOMNAME WII+5   7  C07
ATOMNAME WII+5   8  H01
ATOMNAME WII+5   9  H02
ATOMNAME WII+5   10  H03
ATOMNAME WII+5   11  H04
ATOMNAME WII+5   12  H05
ATOMNAME WII+5   13  C08
ATOMNAME WII+5   14  O01
ATOMNAME WII+5   15  C09
ATOMNAME WII+5   16  C10
ATOMNAME WII+5   17  C11
ATOMNAME WII+5   18  O02
ATOMNAME WII+5   19  N01
ATOMNAME WII+5   20  N02
ATOMNAME WII+5   21  C12
ATOMNAME WII+5   22  C13
ATOMNAME WII+5   23  C14
ATOMNAME WII+5   24  H06
ATOMNAME WII+5   25  H07
ATOMNAME WII+5   26  H08
ATOMNAME WII+5   27  C15
ATOMNAME WII+5   28  H09
ATOMNAME WII+5   29  H10
ATOMNAME WII+5   30  H11
ATOMNAME WII+5   31  C16
ATOMNAME WII+5   32  H12
ATOMNAME WII+5   33  H13
ATOMNAME WII+5   34  H14
ATOMNAME WII+5   35  C17
ATOMNAME WII+5   36  H15
ATOMNAME WII+5   37  H16
ATOMNAME WII+5   38  H17
ATOMNAME WII+5   39  C18
ATOMNAME WII+5   40  H18
ATOMNAME WII+5   41  H19
ATOMNAME WII+5   42  H20
ATOMNAME WII+5   43  C19
ATOMNAME WII+5   44  H21
ATOMNAME WII+5   45  H22
ATOMNAME WII+5   46  H23

CONNECT  WII+5  BR1 ion   
CONNECT  WII+5  C01 ion   0   C02   
CONNECT  WII+5  C02 ion   0   C01   0   C03   0   C07   0   H01   0   H05   
CONNECT  WII+5  C03 ion   0   C02   0   C04   0   H01   0   H02   
CONNECT  WII+5  C04 ion   0   C03   0   C05   0   H01   0   H02   0   H03   
CONNECT  WII+5  C05 ion   0   C04   0   C06   0   H02   0   H03   0   H04   
CONNECT  WII+5  C06 ion   0   C05   0   C07   0   H03   0   H04   0   H05   
CONNECT  WII+5  C07 ion   0   C02   0   C06   0   H04   0   H05   
CONNECT  WII+5  H01 ion   0   C02   0   C03   0   C04   
CONNECT  WII+5  H02 ion   0   C03   0   C04   0   C05   
CONNECT  WII+5  H03 ion   0   C04   0   C05   0   C06   
CONNECT  WII+5  H04 ion   0   C05   0   C06   0   C07   
CONNECT  WII+5  H05 ion   0   C02   0   C06   0   C07   
CONNECT  WII+5  C08 ion   0   O01   
CONNECT  WII+5  O01 ion   0   C08   
CONNECT  WII+5  C09 ion   0   N01   
CONNECT  WII+5  C10 ion   0   N02   
CONNECT  WII+5  C11 ion   0   O02   
CONNECT  WII+5  O02 ion   0   C11   
CONNECT  WII+5  N01 ion   0   C09   0   C13   
CONNECT  WII+5  N02 ion   0   C10   0   C12   
CONNECT  WII+5  C12 ion   0   N02   0   C17   0   H15   0   H16   0   H17   0   C18   0   H18   0   H19   0   H20   0   C19   0   H21   0   H22   0   H23   
CONNECT  WII+5  C13 ion   0   N01   0   C14   0   H06   0   H07   0   H08   0   C15   0   H09   0   H10   0   H11   0   C16   0   H12   0   H13   0   H14   
CONNECT  WII+5  C14 ion   0   C13   0   H06   0   H07   0   H08   
CONNECT  WII+5  H06 ion   0   C13   0   C14   0   H07   0   H08   
CONNECT  WII+5  H07 ion   0   C13   0   C14   0   H06   0   H08   
CONNECT  WII+5  H08 ion   0   C13   0   C14   0   H06   0   H07   
CONNECT  WII+5  C15 ion   0   C13   0   H09   0   H10   0   H11   
CONNECT  WII+5  H09 ion   0   C13   0   C15   0   H10   0   H11   
CONNECT  WII+5  H10 ion   0   C13   0   C15   0   H09   0   H11   
CONNECT  WII+5  H11 ion   0   C13   0   C15   0   H09   0   H10   
CONNECT  WII+5  C16 ion   0   C13   0   H12   0   H13   0   H14   
CONNECT  WII+5  H12 ion   0   C13   0   C16   0   H13   0   H14   
CONNECT  WII+5  H13 ion   0   C13   0   C16   0   H12   0   H14   
CONNECT  WII+5  H14 ion   0   C13   0   C16   0   H12   0   H13   
CONNECT  WII+5  C17 ion   0   C12   0   H15   0   H16   0   H17   
CONNECT  WII+5  H15 ion   0   C12   0   C17   0   H16   0   H17   
CONNECT  WII+5  H16 ion   0   C12   0   C17   0   H15   0   H17   
CONNECT  WII+5  H17 ion   0   C12   0   C17   0   H15   0   H16   
CONNECT  WII+5  C18 ion   0   C12   0   H18   0   H19   0   H20   
CONNECT  WII+5  H18 ion   0   C12   0   C18   0   H19   0   H20   
CONNECT  WII+5  H19 ion   0   C12   0   C18   0   H18   0   H20   
CONNECT  WII+5  H20 ion   0   C12   0   C18   0   H18   0   H19   
CONNECT  WII+5  C19 ion   0   C12   0   H21   0   H22   0   H23   
CONNECT  WII+5  H21 ion   0   C12   0   C19   0   H22   0   H23   
CONNECT  WII+5  H22 ion   0   C12   0   C19   0   H21   0   H23   
CONNECT  WII+5  H23 ion   0   C12   0   C19   0   H21   0   H22   


CHARGE   WII+5 BR1  -0.746610
CHARGE   WII+5 C01  -1.809688
CHARGE   WII+5 C02  1.034025
CHARGE   WII+5 C03  -0.373110
CHARGE   WII+5 C04  -0.204353
CHARGE   WII+5 C05  -0.002636
CHARGE   WII+5 C06  -0.207415
CHARGE   WII+5 C07  -0.366279
CHARGE   WII+5 H01  0.275008
CHARGE   WII+5 H02  0.187741
CHARGE   WII+5 H03  0.148353
CHARGE   WII+5 H04  0.188361
CHARGE   WII+5 H05  0.271753
CHARGE   WII+5 C08  -0.731585
CHARGE   WII+5 O01  0.072381
CHARGE   WII+5 C09  -1.494941
CHARGE   WII+5 C10  -1.498387
CHARGE   WII+5 C11  -0.727602
CHARGE   WII+5 O02  0.070342
CHARGE   WII+5 N01  0.734965
CHARGE   WII+5 N02  0.748416
CHARGE   WII+5 C12  0.449427
CHARGE   WII+5 C13  0.492873
CHARGE   WII+5 C14  -0.531008
CHARGE   WII+5 H06  0.166964
CHARGE   WII+5 H07  0.145541
CHARGE   WII+5 H08  0.122004
CHARGE   WII+5 C15  -0.467070
CHARGE   WII+5 H09  0.110838
CHARGE   WII+5 H10  0.153947
CHARGE   WII+5 H11  0.130220
CHARGE   WII+5 C16  -0.565070
CHARGE   WII+5 H12  0.143982
CHARGE   WII+5 H13  0.189651
CHARGE   WII+5 H14  0.129975
CHARGE   WII+5 C17  -0.536275
CHARGE   WII+5 H15  0.172316
CHARGE   WII+5 H16  0.142399
CHARGE   WII+5 H17  0.131306
CHARGE   WII+5 C18  -0.502994
CHARGE   WII+5 H18  0.117266
CHARGE   WII+5 H19  0.170080
CHARGE   WII+5 H20  0.145591
CHARGE   WII+5 C19  -0.522223
CHARGE   WII+5 H21  0.181005
CHARGE   WII+5 H22  0.124729
CHARGE   WII+5 H23  0.135788
RADIUS   WII+5  BR1   1.85
RADIUS   WII+5  C01   1.70
RADIUS   WII+5  C02   1.70
RADIUS   WII+5  C03   1.70
RADIUS   WII+5  C04   1.70
RADIUS   WII+5  C05   1.70
RADIUS   WII+5  C06   1.70
RADIUS   WII+5  C07   1.70
RADIUS   WII+5  H01   1.20
RADIUS   WII+5  H02   1.20
RADIUS   WII+5  H03   1.20
RADIUS   WII+5  H04   1.20
RADIUS   WII+5  H05   1.20
RADIUS   WII+5  C08   1.70
RADIUS   WII+5  O01   1.52
RADIUS   WII+5  C09   1.70
RADIUS   WII+5  C10   1.70
RADIUS   WII+5  C11   1.70
RADIUS   WII+5  O02   1.52
RADIUS   WII+5  N01   1.55
RADIUS   WII+5  N02   1.55
RADIUS   WII+5  C12   1.70
RADIUS   WII+5  C13   1.70
RADIUS   WII+5  C14   1.70
RADIUS   WII+5  H06   1.20
RADIUS   WII+5  H07   1.20
RADIUS   WII+5  H08   1.20
RADIUS   WII+5  C15   1.70
RADIUS   WII+5  H09   1.20
RADIUS   WII+5  H10   1.20
RADIUS   WII+5  H11   1.20
RADIUS   WII+5  C16   1.70
RADIUS   WII+5  H12   1.20
RADIUS   WII+5  H13   1.20
RADIUS   WII+5  H14   1.20
RADIUS   WII+5  C17   1.70
RADIUS   WII+5  H15   1.20
RADIUS   WII+5  H16   1.20
RADIUS   WII+5  H17   1.20
RADIUS   WII+5  C18   1.70
RADIUS   WII+5  H18   1.20
RADIUS   WII+5  H19   1.20
RADIUS   WII+5  H20   1.20
RADIUS   WII+5  C19   1.70
RADIUS   WII+5  H21   1.20
RADIUS   WII+5  H22   1.20
RADIUS   WII+5  H23   1.20
