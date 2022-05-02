
## Atom hydrophobicity contribution values from Ghose J. Phys. Chem. A 1998, 102, 3762-3772
#
hydrophValuesGhose1998 = {

###################
##     C in      ##
###################

#CH4
1: [1,'[CH4]',-1.5603],

#CH3R
2: [1,'[CH3][#6]',-1.5603],

#CH2R2
3: [2,'[CH2]([#6])[#6]',-1.0120],

#CHR3
4: [3,'[CH]([#6])([#6])[#6]',-0.6681],

#CR4
5: [4,'[C]([#6])([#6])([#6])[#6]',-0.3698],

#CH3X
6: [5,'[CH3X4][N,O,P,S,F,Cl,Br,I]',-1.7880],

#CH2RX
7: [6,'[CH2X4]([#6])[O,N,S,P,Se,F,Cl,Br,I,At]',-1.2486],

#CH2X2
8: [7,'[CH2X4]([N,O,P,S,F,Cl,Br,I])[N,O,P,S,F,Cl,Br,I]',-1.0305],

#CHR2X
9: [8,'[CHX4]([#6])([#6])[O,N,S,P,Se,F,Cl,Br,I,At]',-0.6805],

#CHRX2
10: [9,'[CHX4]([#6])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',-0.3858],

#CHX3
11: [10,'[CHX4]([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.7555],

#CR3X
12: [11,'[CX4]([#6])([#6])([#6])[O,N,S,P,Se,F,Cl,Br,I,At]',-0.2849],

#CR2X2
13: [12,'[CX4]([#6])([#6])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.0200],

#CRX3
14: [13,'[CX4]([#6])([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.7894], 

#CX4
15: [14,'[CX4]([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',1.6422],

#=CH2
16: [15,'[CH2]=[#6]',-0.7866],

#=CHR
17: [16,'[CH1](=[#6])[#6]',-0.3962],

#=CHR
18: [17,'[CH0](=[#6])([#6])[#6]',0.0383],

#=CHX
19: [18,'[CH1](=[#6])[O,N,S,P,Se,F,Cl,Br,I,At]',-0.8051],

#=CRX
20: [19,'[CH0](=[#6])([#6])[O,N,S,P,Se,F,Cl,Br,I,At]',-0.2129],

#=CX2
21: [20,'[CH0](=[#6])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.2432],

#triple-bond-CH
22: [21,'[CH1]#[#6]',0.4697], 

#triple-bond-CR
23: [22,'[CH0](#[#6])[#6]',0.2952], 

#R=C=R
24: [22,'[C](=[#6])=[#6]',0.2952], 

#triple-bond-CH
25: [23,'[CH0](#[#6])[O,N,S,P,Se,F,Cl,Br,I,At]',"undefined"], 

#:R--CH--R aromatic
26: [24,'[cH]',-0.3251], 

#:R--CR--R aromatic
27: [25,'[c](:c)(:c)[#6]',0.1492], 

#:R--CX--R aromatic
28: [26,'[c](:c)(:c)[O,N,S,P,Se,F,Cl,Br,I,At]',0.1539], 

#:R--CH--X aromatic
29: [27,'[cH](:c):[O,N,S,P,Se,F,Cl,Br,I,At]',0.0005], 

#:R--CR--X aromatic
30: [28,'[c]([#6])(:c):[O,N,S,P,Se,F,Cl,Br,I,At]',0.2361], 

#:R--CX--X aromatic
31: [29,'[c]([O,N,S,P,Se,F,Cl,Br,I,At])(:c):[O,N,S,P,Se,F,Cl,Br,I,At]',0.3514], 

#:X--CH--X aromatic
32: [30,'[cH](:[O,N,S,P,Se,F,Cl,Br,I,At]):[O,N,S,P,Se,F,Cl,Br,I,At]',0.1814], 

#:X--CR--X
33: [31,'[c]([#6])(:[O,N,S,P,Se,F,Cl,Br,I,At]):[O,N,S,P,Se,F,Cl,Br,I,At]',0.0901], 

#:X--CX--X
34: [32,'[c]([O,N,S,P,Se,F,Cl,Br,I,At])(:[O,N,S,P,Se,F,Cl,Br,I,At]):[O,N,S,P,Se,F,Cl,Br,I,At]',0.5142],

#:R--CH···X, C···X bond order in pyrrole or furan may be considered as 1
35: [33,'[cr5]([#1])(:[#6]):[O,N,S,P,Se,F,Cl,Br,I,At]',-0.3723], 

#:R--CR···X
36: [34,'[cr5]([#6])(:[#6]):[O,N,S,P,Se,F,Cl,Br,I,At]',0.2813], 

#:R--CX···X 
37: [35,'[cr5]([O,N,S,P,Se,F,Cl,Br,I,At])(:[#6]):[O,N,S,P,Se,F,Cl,Br,I,At]',0.1191], 

#:Al-CH=X, Al represent aliphatic groups
38: [36,'[C]([#1])(=[O,N,S,P,Se,F,Cl,Br,I,At])[C,#1]',-0.1320],

#:Ar-CH=X, Ar represent aromatic groups
39: [37,'[CH1](=[O,N,S,P,Se,F,Cl,Br,I,At])c',-0.0244],

#:Al-C(=X)-Al 
40: [38,'[CH0](=[O,N,S,P,Se,F,Cl,Br,I,At])(C)C',-0.2405],

#:Ar-C(=X)-R
41: [39,'[CH0](=[O,N,S,P,Se,F,Cl,Br,I,At])(c)*',-0.0909],

#:R-C(=X)-X
42: [40,'[CH0](=[O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])*',-0.1002],

#:R-C-triple-bond-X
43: [40,'[CH0](#[O,N,S,P,Se,F,Cl,Br,I,At])*',-0.1002],

#:X=C=X 
44: [40,'[CH0](=[O,N,S,P,Se,F,Cl,Br,I,At])=[O,N,S,P,Se,F,Cl,Br,I,At]',-0.1002],

#:X-C(=X)-X 
45: [41,'[C](=[O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.4182],

#:X--CH···X 
46: [42,'[cr5]([#1])(:[O,N,S,P,Se,F,Cl,Br,I,At]):[O,N,S,P,Se,F,Cl,Br,I,At]',-0.2147],

#:X--CR···X 
47: [43,'[cr5](*)(:[O,N,S,P,Se,F,Cl,Br,I,At]):[O,N,S,P,Se,F,Cl,Br,I,At]',-0.0009],

#:X--CX···X 
48: [44,'[cr5]([O,N,S,P,Se,F,Cl,Br,I,At])(:[O,N,S,P,Se,F,Cl,Br,I,At]):[O,N,S,P,Se,F,Cl,Br,I,At]',0.1388],

#Atom type 45 is unused

###################
## H attached to ##
###################
## The subscript represents hybridization and the superscript its formal oxidation number. The formal oxidation number of a carbon atom equals the sum of the formal bond orders with electronegative atoms; the C--N bond order in pyridine may be considered as 2 while we have one such bond and 1.5 when we have two such bonds; the CâââX bond order in pyrrole or furan may be considered as 1.

# https://www.rdkit.org/docs/RDKit_Book.html
# SMARTS Support and Extensions
# ^0 matches S hybridized atoms
# ^1 matches SP hybridized atoms
# ^2 matches SP2 hybridized atoms
# ^3 matches SP3 hybridized atoms
# ^4 matches SP3D hybridized atoms
# ^5 matches SP3D2 hybridized atoms

#:C^0_sp3 having no X attached to next C 
#49: [46,'[#1][#6;^3][$([#6]~[!O;!N;!S;!P;!Se;!F;!Cl;!Br;!I;!At])]',0.7341],
49: [46,'[#1][#6;^3][#6!$([#6]~N)&!$([#6]~O)&!$([#6]~S)&!$([#6]~P)&!$([#6]~F)&!$([#6]~Cl)&!$([#6]~Br)&!$([#6]~I)]',0.7341],


#:C^1_sp3 
50: [47,'[#1][#6;^3](-[O,N,S,P,Se,F,Cl,Br,I,At])[#6!$([#6]~N)&!$([#6]~O)&!$([#6]~S)&!$([#6]~P)&!$([#6]~F)&!$([#6]~Cl)&!$([#6]~Br)&!$([#6]~I)]',0.6301],

#:C^0_sp2 
#51: [47,'[#1][$([#6;^2]~[!O,!N,!S,!P,!Se,!F,!Cl,!Br,!I,!At])]',0.6301],
51: [47,'[#1][#6;^2]([!O;!N;!S;!P;!Se;!F;!Cl;!Br;!I;!At])[!O;!N;!S;!P;!Se;!F;!Cl;!Br;!I;!At]',0.6301],

#48 :C^2_sp3, C^1_sp2, C^0_sp 0.5180
52: [48,'[#1][$([#6;^3]=[O,N,S,P,Se,F,Cl,Br,I,At])]',0.5180],
53: [48,'[#1][#6;^3]([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.5180],
54: [48,'[#1][$([#6;^2]-[O,N,S,P,Se,F,Cl,Br,I,At])]',0.5180],
55: [48,'[#1][$([#6;^1]~[!O,!N,!S,!P,!Se,!F,!Cl,!Br,!I,!At])]',0.5180],

#49 :C^3_sp3, C^2-3_sp2, C^1-3_sp -0.0371
56: [49,'[#1][$([#6;^3]#[O,N,S,P,Se,F,Cl,Br,I,At])]',-0.0371],
57: [49,'[#1][$([#6;^3](=[O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At])]',-0.0371],
58: [49,'[#1][$([#6;^3]([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At])]',-0.0371],
59: [49,'[#1][$([#6;^2]=[O,N,S,P,Se,F,Cl,Br,I,At])]',-0.0371],
60: [49,'[#1][$([#6;^2]([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At])]',-0.0371],
61: [49,'[#1][$([#6;^1][O,N,S,P,Se,F,Cl,Br,I,At])]',-0.0371],

#:heteroatom
62: [50,'[#1][*;!#6]',-0.1036],

# 51 :alpha-Cd, An R-C-alpha may be defined as a C attached through a single bond with -C=X, -C-triple-bond-X, -C--X aromatic bond
63: [51,'[#1][#6]-[$([#6]=[O,N,S,P,Se,F,Cl,Br,I,At])]',0.5234],
64: [51,'[#1][#6]-[$([#6]#[O,N,S,P,Se,F,Cl,Br,I,At])]',0.5234],
65: [51,'[#1][#6]-[$([#6]:[O,N,S,P,Se,F,Cl,Br,I,At])]',0.5234],

#52 :C0_sp3, having 1 X attached to next carbon
66: [52,'[#1][#6;^3][$([#6][O,N,S,P,Se,F,Cl,Br,I,At])]',0.6666],

#53 :C0_sp3, having 2 X attached to next carbon
67: [53,'[#1][#6;^3][$([#6]([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At])]',0.5372],

#54 :C0_sp3, having 3 X attached to next carbon
68: [54,'[#1][#6;^3][$([#6]([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At])]',0.6338],

#55 :C0_sp3, having 4 or more X attached to next carbon
69: [55,'[#1][#6;^3][$([#6]([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])([O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At])]',0.3620],

###################
##     O in      ##
###################

#56 :alcohol
70: [56,'[OX2H]-[!c;!#6X3;!$([CX3]=[OX1])]',-0.3567],

#57 :phenol, enol, carboxyl OH -0.0127
71: [57,'[OX2]([H])[cX3]:[c]',-0.0127],
72: [57,'[OX2]([H])[#6X3]=[#6]',-0.0127],
73: [57,'[OX2]([H])[$([CX3]=[OX1])]',-0.0127],


#58 :=O -0.0233
74: [58,'[OX1]=[#6]',-0.0233],

#59 :Al-O-Al -0.1541, Al represent aliphatic groups
75: [59,'[OX2](C)C',-0.1541],

#60 :Al-O-Ar, Ar2O, :R···O···R, R-O-C=X 0.0324, Ar represent aromatic groups
76: [60,'[OX2](C)c',0.0324],
77: [60,'[OX2](c)c',0.0324],
78: [60,'[OX2]([#6])[#6]',0.0324],
79: [60,'[OX2]([#6])[$([#6]=[O,N,S,P,Se,F,Cl,Br,I,At])]',0.0324],

#61 :--O 1.0520, as in nitro, N-oxides.
80: [61,'[OX1-][#7+]',1.0520],
81: [61,'[OX1]=[#7]',1.0520],

#62 :O- (negatively charged) -0.7941
82: [62,'[OX1-][*;#7+]',-0.7941],

#63 :R-O-O-R 0.4165
83: [63,'[OX2][OX2]',0.4165],

###################
##    Se in      ##
###################

#64 :Any-Se-Any 0.6601
84: [64,'[SeX2](*)*',0.6601],

#65 :=Se, not used

###################
##     N in      ##
###################

#66 :Al-NH2 -0.5427, Al represent aliphatic groups
85: [66,'[NX3]([#1])([#1])A',-0.5427],

#67 :Al2NH -0.3168
86: [67,'[NX3]([#1])(A)A',-0.3168],

#68 :Al3N 0.0132
87: [68,'[NX3](A)(A)A',0.0132],

#69 :Ar-NH2, X-NH2 -0.3883, Ar represent aromatic groups
88: [69,'[NX3]([#1])([#1])a',-0.3883],
89: [69,'[NX3]([#1])([#1])[O,N,S,P,Se,F,Cl,Br,I,At]',-0.3883],

#70 :Ar-NH-Al -0.0389
90: [70,'[NX3]([#1])(A)a',-0.0389],

#71 :Ar-NAl2 0.1087
91: [71,'[NX3](A)(A)a',0.1087],

#72 :RCO-N<, >N-X=X -0.5113
92: [72,'[NX3](*)(*)[$(C(=O)*)]',-0.5113],
93: [72,'[NX3](*)(*)[$([O,N,S,P,Se,F,Cl,Br,I,At]=[O,N,S,P,Se,F,Cl,Br,I,At])]',-0.5113],

#73 :Ar2NH, Ar3N, Ar2N-Al, R···N···R (Pyrrole-type structure) 0.1259
94: [73,'[NX3]([#1])(a)a',0.1259],
95: [73,'[NX3](a)(a)a',0.1259],
96: [73,'[NX3](A)(a)a',0.1259],
97: [73,'[NR5X2](*)*',0.1259],

#74 :R-triple-bond-N, R=N- 0.1349
98: [74,'[NX1]#*',0.1349],
99: [74,'[NX2]=*',0.1349],

#75 :R--N--R, R--N--X (Pyridine type structure) -0.1624
100: [75,'[NX2](*)*',-0.1624],
101: [75,'[NX2]([O,N,S,P,Se,F,Cl,Br,I,At])*',-0.1624],

#76 :Ar-NO2, R--N(--R)--O (Pyridine N-oxide type structure), RO-NO, -2.0585
102: [76,'[NX3](O)(O)a',-2.0585],
103: [76,'[NX3](O)(*)*',-2.0585],
104: [76,'[NX2](O)[$(NO)]',-2.0585],

#77 :Al-NO2 -1.9150
105: [77,'[NX3](O)(O)A',-1.9150],

#78 :Ar-N=X, X-N=X 0.4208
106: [78,'[NX2](=[O,N,S,P,Se,F,Cl,Br,I,At])A',0.4208],
107: [78,'[NX2](=[O,N,S,P,Se,F,Cl,Br,I,At])[O,N,S,P,Se,F,Cl,Br,I,At]',0.4208],

#79 :N+ (positively charged) -1.4439
108: [79,'[#7+]',-1.4439],
109: [79,'[#7X4]',-1.4439],
109: [79,'[#7X3](=*)*',-1.4439],

#80 unused



} # End hydrophValues dict


## Atom-based calculation of LogP and MR using Crippen's approach
# Adapted from RDKit Crippen.txt file
#    Reference:
#      S. A. Wildman and G. M. Crippen *JCICS* _39_ 868-873 (1999)

hydrophValuesCrippen1999={
1:[1,'[CH4]',0.1441],#C1
2:[2,'[CH3]C',0.1441],#C1
3:[3,'[CH2](C)C',0.1441],#C1
4:[4,'[CH](C)(C)C',0],#C2
5:[5,'[C](C)(C)(C)C',0],#C2
6:[6,'[CH3][N,O,P,S,F,Cl,Br,I]',-0.2035],#C3
7:[7,'[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]',-0.2035],#C3
8:[8,'[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]',-0.2051],#C4
9:[9,'[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]',-0.2051],#C4
10:[10,'[C]=[!C;A;!#1]',-0.2783],#C5
11:[11,'[CH2]=C',0.1551],#C6
12:[12,'[CH1](=C)[A;!#1]',0.1551],#C6
13:[13,'[CH0](=C)([A;!#1])[A;!#1]',0.1551],#C6
14:[14,'[C](=C)=C',0.1551],#C6
15:[15,'[CX2]#[A;!#1]',0.0017],#C7
16:[16,'[CH3]c',0.08452],#C8
17:[17,'[CH3]a',-0.1444],#C9
18:[18,'[CH2X4]a',-0.0516],#C10
19:[19,'[CHX4]a',0.1193],#C11
20:[20,'[CH0X4]a',-0.0967],#C12
21:[21,'[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]',-0.5443],#C13
22:[22,'[c][#9]',0],#C14
23:[23,'[c][#17]',0.245],#C15
24:[24,'[c][#35]',0.198],#C16
25:[25,'[c][#53]',0],#C17
26:[26,'[cH]',0.1581],#C18
27:[27,'[c](:a)(:a):a',0.2955],#C19
28:[28,'[c](:a)(:a)-a',0.2713],#C20
29:[29,'[c](:a)(:a)-C',0.136],#C21
30:[30,'[c](:a)(:a)-N',0.4619],#C22
31:[31,'[c](:a)(:a)-O',0.5437],#C23
32:[32,'[c](:a)(:a)-S',0.1893],#C24
33:[33,'[c](:a)(:a)=[C,N,O]',-0.8186],#C25
34:[34,'[C](=C)(a)[A;!#1]',0.264],#C26
35:[35,'[C](=C)(c)a',0.264],#C26
36:[36,'[CH1](=C)a',0.264],#C26
37:[37,'[C]=c',0.264],#C26
38:[38,'[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]',0.2148],#C27
39:[39,'[#6]',0.08129],#CS
40:[40,'[#1][#6,#1]',0.123],#H1
41:[41,'[#1]O[CX4,c]',-0.2677],#H2
42:[42,'[#1]O[!C;!N;!O;!S]',-0.2677],#H2
43:[43,'[#1][!C;!N;!O]',-0.2677],#H2
44:[44,'[#1][#7]',0.2142],#H3
45:[45,'[#1]O[#7]',0.2142],#H3
46:[46,'[#1]OC=[#6,#7,O,S]',0.298],#H4
47:[47,'[#1]O[O,S]',0.298],#H4
48:[48,'[#1]',0.1125],#HS
49:[49,'[NH2+0][A;!#1]',-1.019],#N1
50:[50,'[NH+0]([A;!#1])[A;!#1]',-0.7096],#N2
51:[51,'[NH2+0]a',-1.027],#N3
52:[52,'[NH1+0]([!#1;A,a])a',-0.5188],#N4
53:[53,'[NH+0]=[!#1;A,a]',0.08387],#N5
54:[54,'[N+0](=[!#1;A,a])[!#1;A,a]',0.1836],#N6
55:[55,'[N+0]([A;!#1])([A;!#1])[A;!#1]',-0.3187],#N7
56:[56,'[N+0](a)([!#1;A,a])[A;!#1]',-0.4458],#N8
57:[57,'[N+0](a)(a)a',-0.4458],#N8
58:[58,'[N+0]#[A;!#1]',0.01508],#N9
59:[59,'[NH3,NH2,NH;+,+2,+3]',-1.95],#N10
60:[60,'[n+0]',-0.3239],#N11
61:[61,'[n;+,+2,+3]',-1.119],#N12
62:[62,'[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]',-0.3396],#N13
63:[63,'[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]',-0.3396],#N13
64:[64,'[NH0;+,+2,+3](=[#6])=[#7]',-0.3396],#N13
65:[65,'[N;+,+2,+3]#[A;!#1]',0.2887],#N14
66:[66,'[N;-,-2,-3]',0.2887],#N14
67:[67,'[N;+,+2,+3](=[N;-,-2,-3])=N',0.2887],#N14
68:[68,'[#7]',-0.4806],#NS
69:[69,'[o]',0.1552],#O1
70:[70,'[OH,OH2]',-0.2893],#O2
71:[71,'[O]([A;!#1])[A;!#1]',-0.0684],#O3
72:[72,'[O](a)[!#1;A,a]',-0.4195],#O4
73:[73,'[O]=[#7,#8]',0.0335],#O5
74:[74,'[OX1;-,-2,-3][#7]',0.0335],#O5
75:[75,'[OX1;-,-2,-2][#16]',-0.3339],#O6
76:[76,'[O;-0]=[#16;-0]',-0.3339],#O6
77:[77,'[O-]C(=O)',-1.326],#O12
78:[78,'[OX1;-,-2,-3][!#1;!N;!S]',-1.189],#O7
79:[79,'[O]=c',0.1788],#O8
80:[80,'[O]=[CH]C',-0.1526],#O9
81:[81,'[O]=C(C)([A;!#1])',-0.1526],#O9
82:[82,'[O]=[CH][N,O]',-0.1526],#O9
83:[83,'[O]=[CH2]',-0.1526],#O9
84:[84,'[O]=[CX2]=O',-0.1526],#O9
85:[85,'[O]=[CH]c',0.1129],#O10
86:[86,'[O]=C([C,c])[a;!#1]',0.1129],#O10
87:[87,'[O]=C(c)[A;!#1]',0.1129],#O10
88:[88,'[O]=C([!#1;!#6])[!#1;!#6]',0.4833],#O11
89:[89,'[#8]',-0.1188],#OS
90:[90,'[#9-0]',0.4202],#F
91:[91,'[#17-0]',0.6895],#Cl
92:[92,'[#35-0]',0.8456],#Br
93:[93,'[#53-0]',0.8857],#I
94:[94,'[#9,#17,#35,#53;-]',-2.996],#Hal
95:[95,'[#53;+,+2,+3]',-2.996],#Hal
96:[96,'[+;#3,#11,#19,#37,#55]',-2.996],#Hal
97:[97,'[#15]',0.8612],#P
98:[98,'[S;-,-2,-3,-4,+1,+2,+3,+5,+6]',-0.0024],#S2
99:[99,'[S-0]=[N,O,P,S]',-0.0024],#S2
100:[100,'[S;A]',0.6482],#S1
101:[101,'[s;a]',0.6237],#S3
102:[102,'[#3,#11,#19,#37,#55]',-0.3808],#Me1
103:[103,'[#4,#12,#20,#38,#56]',-0.3808],#Me1
104:[104,'[#5,#13,#31,#49,#81]',-0.3808],#Me1
105:[105,'[#14,#32,#50,#82]',-0.3808],#Me1
106:[106,'[#33,#51,#83]',-0.3808],#Me1
107:[107,'[#34,#52,#84]',-0.3808],#Me1
108:[108,'[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]',-0.0025],#Me2
109:[109,'[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]',-0.0025],#Me2
110:[110,'[#72,#73,#74,#75,#76,#77,#78,#79,#80]',-0.0025],#Me2

}
