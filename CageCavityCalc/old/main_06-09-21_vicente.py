import numpy as np
import math
import time
from scipy.spatial import KDTree
import MDAnalysis
import os
from functools import partial
from multiprocess import Pool

import rdkit.Chem.AllChem as rdkit
from rdkit.Geometry import Point3D
from rdkit import Chem

# there are a lot of warning, since not everythin is specify in pdb, let's ignore them
import warnings
warnings.filterwarnings("ignore")

# Van der Waals radii in Å taken from http://www.webelements.com/periodicity/van_der_waals_radius/

vdw_radii = {'h': 1.20, 'he': 1.40, 'li': 1.82, 'be': 1.53, 'b': 1.92, 'c': 1.70, 'n': 1.55, 'o': 1.52, 'f': 1.47,
             'ne': 1.54, 'na': 2.27, 'mg': 1.73, 'al': 1.84, 'si': 2.10, 'p': 1.80, 's': 1.80, 'cl': 1.75, 'ar': 1.88,
             'k': 2.75, 'ca': 2.31, 'ni': 1.63, 'cu': 1.40, 'zn': 1.39, 'ga': 1.87, 'ge': 2.11, 'as': 1.85, 'se': 1.90,
             'br': 1.85, 'kr': 2.02, 'rb': 3.03, 'sr': 2.49, 'pd': 1.63, 'ag': 1.72, 'cd': 1.58, 'in': 1.93, 'sn': 2.17,
             'sb': 2.06, 'te': 2.06, 'i': 1.98, 'xe': 2.16, 'cs': 3.43, 'ba': 2.49, 'pt': 1.75, 'au': 1.66, 'fe':1.40}

atoms_lib = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'sc',
             'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y',
             'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce',
             'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os',
             'ir', 'pt', 'au',  'hg', 'tl', 'pb', 'bi','po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu',
             'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn',
             'nh', 'fl', 'mc', 'lv', 'ts', 'og']



## Atom hydrophobicity contribution values from Ghose J. Phys. Chem. A 1998, 102, 3762-3772
#
hydrophValues = {

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
108: [79,'[N+]',-1.4439],

#80 unused



} # End hydrophValues dict



def cavity_volume(positions,radius=1.2, volume_grid_size=0.2):
    # Find min, max for the grid box, and create x,y,z of grid
    x_space = np.arange(np.min(positions.T[0])-radius, np.max(positions.T[0])+radius, volume_grid_size)
    y_space = np.arange(np.min(positions.T[1])-radius, np.max(positions.T[1])+radius, volume_grid_size)
    z_space = np.arange(np.min(positions.T[2])-radius, np.max(positions.T[2])+radius, volume_grid_size)
    # convert x,y,z axis to 3d vector
    grid_points = np.vstack(np.meshgrid(x_space,y_space,z_space)).reshape(3,-1).T
    # calcualte distance between the points
    dist_matrix = MDAnalysis.lib.distances.distance_array(positions, grid_points)
    # if grid point within radius of position then we add volume of the grid
    return np.sum(np.sum(dist_matrix<radius,axis=0)>0)*(volume_grid_size**3)

def cavity(frame_index, syst, intput, output, grid_spacing = 1.0, distance_threshold_for_90_deg_angle = 7, save_only_surface = True,calculate_bfactor = True, compute_aromatic_contacts = False, compute_atom_contacts = True, distThreshold_atom_contacts = 5.0, number_of_dummies=500, dummy_atom_diameter=vdw_radii['h'], pore_radius_limit=10):

    syst.universe.trajectory[frame_index]
    if frame_index>0:
        print("Progress: {:d}/{:d}".format(frame_index, len(syst.universe.trajectory)), end = '\r')

    ########## cavity calculation
    cagePDB = intput
    # we need mol files as pdb file was not behaving OK with the aromatic
    cageMOL = cagePDB.replace(".pdb", ".mol2")
    #print(cagePDB)
    cagePDBout1 = cagePDB.replace(".pdb", ".cavity.pdb")

    calculate_windows = True
    #threads_KDThree = 4

    ######################################################################
    ###   Calculate cavity
    ######################################################################
    #print("\n--- Calculation of the cavity ---")

    start_time = time.time()

    ## Using the PDB file as the MOL file for this example is not working, some issues as it was generated from the CIF file
    rdkit_cage = rdkit.MolFromMol2File(cageMOL, removeHs=False)
    numberOfAtoms = rdkit_cage.GetNumAtoms()
    
        
    listGlobalOut = list()
    listGlobalOut2 = list() 
    for i in hydrophValues:
        #print("Dict list number " + str(i) + " " + hydrophValues[i][1])
        match = rdkit_cage.GetSubstructMatches(Chem.MolFromSmarts(hydrophValues[i][1]),False,False)
        #print(match)
        listOut = list()
        for x in match:
            #print(x[0]+1)
            listOut.append(x[0]+1)
        
        listGlobalOut2.append(list(set(listOut)))
        if listOut:
            #print(listOut)
            #print(list(set(listOut)))
            listGlobalOut.extend(list(set(listOut)))

    #print(len(listGlobalOut2))
    #print(listGlobalOut2)
    #print(list(set(listGlobalOut)))
    #print("++++++++++++++++++++++\n\n")

    numberOfAtoms = rdkit_cage.GetNumAtoms()

    atomTypesList = []
    for i in range(0,numberOfAtoms):
        atomTypesList.append([])
    #print(atomTypesList)

    atomTypesListHydrophValues = []
    for i in range(0,numberOfAtoms):
        atomTypesListHydrophValues.append([])
    #print(atomTypesList)

    atomTypesMeanListHydrophValues = []

    for i in range(0, len(listGlobalOut2)):
        #print("Atom type " + str(i))
        if listGlobalOut2[i]:
            #print(i)
            for j in listGlobalOut2[i]:
                #print(j)            
                atomTypesList[j-1].append(i+1)
                #print(atomTypesList[j-1])


    for i in range(0,numberOfAtoms):
        atomSymbol = rdkit_cage.GetAtomWithIdx(i).GetSymbol()
        valuesList = []
        if len(atomTypesList[i])>0:
            valuesList = []
            for j in atomTypesList[i]:
                valuesList.append(hydrophValues[j][2])
        
        atomTypesListHydrophValues.append(valuesList)
        meanValuestList = np.mean(valuesList)
        atomTypesMeanListHydrophValues.append(meanValuestList)
        print(atomSymbol, i+1, atomTypesList[i], valuesList, meanValuestList)


    #print(atomTypesMeanListHydrophValues)
    
    
    '''
    rdkit_cage = rdkit.MolFromPDBFile(cagePDB,sanitize=False, removeHs=False, proximityBonding=False)
    
    # Calculate the centeor of mass
    atoms = [atom for atom in rdkit_cage.GetAtoms()]
    numberOfAtoms = rdkit_cage.GetNumAtoms()
    xyzPositionArray = np.array([list(rdkit_cage.GetConformer().GetAtomPosition(i)) for i in range(numberOfAtoms)])
    center_of_mass = np.array(
        sum(atoms[i].GetMass() * xyzPositionArray[i] for i in range(numberOfAtoms))) / Descriptors.MolWt(rdkit_cage)
    print("center_of_mass= ", center_of_mass)

    # Get the xyx atom coordinates
    xyzAtoms = np.array(
        [rdkit_cage.GetConformer().GetAtomPosition(cage_atom.GetIdx()) for cage_atom in rdkit_cage.GetAtoms()],
        dtype=np.float64)
    '''
    center_of_mass = syst.atoms.center_of_mass()
    #print("center_of_mass= ", center_of_mass)
    xyzAtoms = syst.atoms.positions

    # Calculate neighbors with KDTree. Leaf_size affects the speed of a query and the memory required to store the constructed tree.
    # The amount of memory needed to store the tree scales as approximately n_samples / leaf_size.
    kdtxyzAtoms = KDTree(xyzAtoms, leafsize=20)
    # Calculate the set of atoms within the euclidean distanc upper bound
    result = kdtxyzAtoms.query(xyzAtoms[30], k=None, p=2, distance_upper_bound=3.0)

    # We use the center of mass of the molecule as the pore center of mass and the pore radius as the 95% of the close contact to the center_of_mass
    pore_center_of_mass = center_of_mass
    distancesFromCOM = kdtxyzAtoms.query(center_of_mass, k=None, p=2)
    pore_radius = distancesFromCOM[0][0] * 0.8
    maxDistanceFromCOM = (distancesFromCOM[0][-1])
    pore_diameter = pore_radius * 2

    class GridPoint:
        def __init__(self, i, j, k, points, center, grid_spacing):
            self.i = i
            self.j = j
            self.k = k

            self.i_from_center = i - points
            self.j_from_center = j - points
            self.k_from_center = k - points

            self.x = self.i_from_center * grid_spacing + center[0]
            self.y = self.j_from_center * grid_spacing + center[1]
            self.z = self.k_from_center * grid_spacing + center[2]

            self.pos = [self.x, self.y, self.z]

            self.d_from_center = np.linalg.norm(np.array(self.pos) - np.array(center))
            self.inside_cavity = 0
            self.overlapping_with_cage = 0
            self.number_of_neighbors = 0
            self.is_window = 0
            self.vector_angle = 0
            self.neighbors = []

    class CageGrid:
        def __init__(self, center, radius, delta, grid_spacing):
            self.center = center
            self.size = radius + delta
            self.grid_spacing = grid_spacing
            self.points = math.ceil(self.size / self.grid_spacing)
            self.n = 2 * self.points + 1
            self.grid = [GridPoint(i, j, k, self.points, self.center, self.grid_spacing) for k in range(self.n) for j in
                         range(self.n) for i in range(self.n)]
            self.grid = np.array(self.grid)
            self.gridPosList = []
            for i in self.grid:
                self.gridPosList.append(i.pos)

    box_size = maxDistanceFromCOM
    #print("box_size=", box_size)
    calculatedGird = CageGrid(pore_center_of_mass, box_size, delta=0, grid_spacing=grid_spacing)

    # -#-#-#-# KDTree algorithm to speed up the method, needs to be implemented with the code below

    # Create a KDThree of the dummy atom gird
    calculatedGirdKDTree = KDTree(calculatedGird.gridPosList, leafsize=20)
    # Calculate contacts of dummy atoms with a distance < 1.1*grid_spacing, rages from 1 to 7 (as 1 contact is always, the atom itself)

    '''
    for dummy_atom in calculatedGird.grid:
        xyzDummySet = calculatedGirdKDTree.query(dummy_atom.pos, k=None, p=2, distance_upper_bound=1.1*grid_spacing)
        dummy_atom.neighbors = xyzDummySet[1] #does not work as the index does not correspond
    '''

    ## IMPROVED ALGORITHM USING KDThree
    atom_type_list = []
    coords_dict = {}
    #vdwR_dict = {}
    atom_idx_dict = {}

    for cage_atom in syst.atoms:
        if cage_atom.name.lower()[:2] in atoms_lib:
            atom_type = cage_atom.name[:2].lower()
        else:
            atom_type = cage_atom.name[0].lower()

        atom_idx = cage_atom.index
        #pos = rdkit_cage.GetConformer().GetAtomPosition(atom_idx)
        pos = cage_atom.position


        if atom_type not in atom_type_list:
            atom_type_list.append(atom_type)
        '''
            #vdwR = rdkit.GetPeriodicTable().GetRvdw(cage_atom.GetAtomicNum())
            vdwR = vdw_radii[]
            vdwR_dict.setdefault(atom_type, []).append(vdwR)
        '''
        coords_dict.setdefault(atom_type, []).append(list(pos))
        atom_idx_dict.setdefault(atom_type, []).append(atom_idx)

    # print(atom_type_list)
    # print(vdwR_dict)
    # print(coords_dict)
    # print(atom_idx_dict)

    # Create a KDTree for each atom type and store in a dict
    KDTree_dict = {}
    for atom_type in atom_type_list:
        # print(coords_dict[atom_type])
        # print(atom_type)
        KDTree_dict[atom_type] = KDTree(coords_dict[atom_type], leafsize=20)

    # Radius of the dummy H atoms in the cavity
    inside_cavity_grid = []
    #vdwRdummy = vdw_radii['h']# rdkit.GetPeriodicTable().GetRvdw(1)
    vdwRdummy = dummy_atom_diameter
    for a, i in enumerate(calculatedGird.grid):
        dist1 = np.linalg.norm(np.array(pore_center_of_mass) - np.array(i.pos))
        if dist1 > 0:
            vect1Norm = (np.array(pore_center_of_mass) - np.array(i.pos)) / dist1
        #if (dist1 < pore_radius and pore_radius > 10):
        if (dist1 < pore_radius_limit):
            i.inside_cavity = 1
            inside_cavity_grid.append(a)
        if i.inside_cavity == 0:
            summAngles = []
            distancesAngles = []
            for atom_type in atom_type_list:
                xyzAtomsSet2 = KDTree_dict[atom_type].query(i.pos, k=None, p=2,
                                                            distance_upper_bound=distance_threshold_for_90_deg_angle)
                if xyzAtomsSet2[1]:
                    #vdwR = vdwR_dict[atom_type][0]
                    vdwR = vdw_radii[atom_type]
                    distThreshold = vdwR + vdwRdummy
                    if xyzAtomsSet2[0][0] < distThreshold:
                        # print("Dummy atom overlapping with cage")
                        i.overlapping_with_cage = 1
                        # break

                    # check that distances are the same by KDTree query an position based calculation
                    """"
                    for z in range (0,len(xyzAtomsSet2[1])):
                        atom_pos = coords_dict[atom_type][xyzAtomsSet2[1][z]]
                        dist2 = np.linalg.norm(np.array(atom_pos)-np.array(i.pos))
                        print(xyzAtomsSet2[0][z], dist2)
                    """

                    for atom_pos_index in xyzAtomsSet2[1]:
                        atom_pos = coords_dict[atom_type][atom_pos_index]
                        dist2 = np.linalg.norm(np.array(atom_pos) - np.array(i.pos))
                        vect2Norm = (np.array(atom_pos) - np.array(i.pos)) / dist2
                        angle = np.arccos(np.dot(vect1Norm, vect2Norm))
                        summAngles.append(angle)
                        distancesAngles.append(1 / dist2)

            if (summAngles):
                # averageSummAngles = sum(summAngles) / len(summAngles)
                averageSummAngles = np.average(summAngles, axis=None, weights=distancesAngles)
                averageSummAngles_deg = np.degrees(averageSummAngles)
                i.vector_angle = averageSummAngles_deg
                if i.overlapping_with_cage == 0:
                    if averageSummAngles_deg > 90:
                        # print("Added dummy atom to cavity", np.degrees(averageSummAngles), i.overlapping_with_cage)
                        i.inside_cavity = 1
                        inside_cavity_grid.append(a)
            """
            for j, cage_atom in enumerate(rdkit_cage.GetAtoms()):
            #for cage_atom in rdkit_cage.GetAtoms():
                if len(xyzAtomsSet[1]):
                    if (cage_atom.GetIdx() in xyzAtomsSet[1]):
                        vdwR = rdkit.GetPeriodicTable().GetRvdw(cage_atom.GetAtomicNum())
                        distThreshold = vdwR + vdwRdummy
                        xyzAtomsSet2 = kdtxyzAtoms.query(i.pos, k = None, p = 2, distance_upper_bound = distThreshold)
                        if xyzAtomsSet2[1] == []:
                            pos1 = rdkit_cage.GetConformer().GetAtomPosition(cage_atom.GetIdx())                
                            dist2 = np.linalg.norm(np.array(pos1)-np.array(i.pos))
                            vect2Norm = (np.array(pos1)-np.array(i.pos)) / dist2
                            angle = np.arccos(np.dot(vect1Norm, vect2Norm))

                            summAngles.append(angle)
                            distancesAngles.append(1/dist2)

                            #if np.degrees(angle) < 150:
                            #    summAngles.append(angle)

                            #if np.degrees(angle) > 140:
                            #    i.inside_cavity = 1
                            #    break
            """

    # Create a KDThree of the dummy atom girds to  remove isolated dummy atoms with a distance < 2*vdwRdummy
    calculatedGirdContacts = []
    for i in calculatedGird.grid[inside_cavity_grid]:
            calculatedGirdContacts.append(i.pos)
    calculatedGirdContactsKDTree = KDTree(calculatedGirdContacts, leafsize=20)

    # Calculate contacts of dummy atoms with a distance < 1.1*grid_spacing, rages from 1 to 7 (as 1 contact is alwas, the atom itself)
    for i, dummy_atom in enumerate(calculatedGird.grid[inside_cavity_grid]):
            xyzDummySet = calculatedGirdContactsKDTree.query(dummy_atom.pos, k=None, p=2,
                                                             distance_upper_bound=1.1 * grid_spacing)
            dummy_atom.number_of_neighbors = len(xyzDummySet[1])

    # Create a KDThree of the dummy atom girds to that overlap with the cage

    for i in calculatedGird.grid:
        for atom_type in atom_type_list:
            #vdwR = vdwR_dict[atom_type][0]
            vdwR = vdw_radii[atom_type]
            distThreshold = vdwR + vdwRdummy
            xyzAtomsSet2 = KDTree_dict[atom_type].query(i.pos, k=None, p=2, distance_upper_bound=distThreshold)
            if xyzAtomsSet2[0]:
                i.overlapping_with_cage = 1
            else:
                i.overlapping_with_cage = 0

    overlappingCalculatedGirdContacts = []
    overlappingCalculatedGirdContacts_angles = []
    for i in calculatedGird.grid:
        if i.overlapping_with_cage == 1:
            # if i.vector_angle > 90:
            overlappingCalculatedGirdContacts.append(i.pos)
            overlappingCalculatedGirdContacts_angles.append(i.vector_angle)
    overlappingCalculatedGirdContactsKDTree = KDTree(overlappingCalculatedGirdContacts, leafsize=20)

    ## IMPROVED ALGORITHM USING KDThree
    # Remove dummy atoms overapping by the cage

    # distance111 = KDTree_dict['Pd'].query([0,0,0], k = None, p = 2, distance_upper_bound = 50)
    # print(distance111)
    """
    for i, dummy_atom in enumerate(calculatedGird.grid):
        if dummy_atom.inside_cavity == 1:
            for atom_type in atom_type_list:
                vdwR = vdwR_dict[atom_type][0]
                distThreshold = vdwR + vdwRdummy
                dist1 = KDTree_dict[atom_type].query(dummy_atom.pos, k = 1, p = 2)
                if dist1[0] < 1.01*distThreshold:
                    dummy_atom.inside_cavity = 0
    """
    """
    for i, dummy_atom in enumerate(calculatedGird.grid):
        if dummy_atom.inside_cavity == 1:
            for cage_atom in rdkit_cage.GetAtoms():
                pos1 = rdkit_cage.GetConformer().GetAtomPosition(cage_atom.GetIdx())
                vdwR = rdkit.GetPeriodicTable().GetRvdw(cage_atom.GetAtomicNum())
                distThreshold = vdwR + vdwRdummy
                dist1 = np.linalg.norm(np.array(pos1)-np.array(dummy_atom.pos))
                if dist1 < 1.01*distThreshold:
                    dummy_atom.inside_cavity = 0
    """

    # Create an empty editable molecule
    if (number_of_dummies == 0):
        number_of_dummies = len(calculatedGird.grid)

    dummy_universe = MDAnalysis.Universe.empty(number_of_dummies, n_residues=number_of_dummies,
                                               atom_resindex=[0] * number_of_dummies,
                                               residue_segindex=[0] * number_of_dummies, trajectory=True)
    dummy_universe.add_TopologyAttr('tempfactors')
    # we put dummies far away from cage:
    dummy_universe.atoms.positions = number_of_dummies * np.array([0,0,0])
    dummy_universe.add_TopologyAttr('name', ['D'] * number_of_dummies)
    dummy_universe.add_TopologyAttr('resname', ['D'] * number_of_dummies)
        
    dummies_inside = []

    for i, dummy_atom in enumerate(calculatedGird.grid[inside_cavity_grid]):
        if (dummy_atom.inside_cavity == 1 and dummy_atom.overlapping_with_cage == 0):
            if calculate_windows == True:
                distWindow = overlappingCalculatedGirdContactsKDTree.query(dummy_atom.pos, k=None, p=2,
                                                                           distance_upper_bound=1.1 * grid_spacing)
                # print(distWindow)
                if len(distWindow[0]) == 0:
                    dummy_atom.is_window = 1
                if len(distWindow[0]) > 0:
                    if overlappingCalculatedGirdContacts_angles[
                        distWindow[1][0]] < 100:  ## 100 deg is necessary to include the atoms close to the surface
                        dummy_atom.is_window = 1
                        # print("dummy atom not in a window")

            # Remove isolated dummy atoms with a distance < 1.1*grid_spacing
            #if (dummy_atom.number_of_neighbors > 0):
            if (dummy_atom.number_of_neighbors > 1 and dummy_atom.number_of_neighbors < 7 and save_only_surface == True):
                if calculate_bfactor == True:
                    contacts = []
                    for atom_type in atom_type_list:
                        if atom_type:
                        #if atom_type == "c":
                        #if atom_type != "h":
                            dist = KDTree_dict[atom_type].query(dummy_atom.pos, k=None, p=2,
                                                                distance_upper_bound=distThreshold_atom_contacts)
                            
                            #print(dist)

                            if dist[0]:
                                for k in range(0, len(dist[0])):
                                    if compute_atom_contacts == True:
                                        # print("compute_atom_contacts")
                                        #contacts.append(1)
                                        #contacts.append(1 / dist[0][k])
                                        
                                        #Fauchtre's distance function exp(-d), Journal of Computer-Aided Molecular Design, 1994, 8, 83-96
                                        #contacts.append(atomTypesMeanListHydrophValues[dist[1][k]]*np.exp(-1*dist[0][k]))
                                        
                                        #Audry's distance function (1/(1+d)), Journal of Computer-Aided Molecular Design, 1994, 8, 83-96
                                        contacts.append(atomTypesMeanListHydrophValues[dist[1][k]]/(1+dist[0][k]))
                                        #print(atomTypesMeanListHydrophValues[dist[1][k]])
                                        
                                    '''
                                    elif compute_aromatic_contacts == True:
                                        isAromatic = rdkit_cage.GetAtomWithIdx(
                                            atom_idx_dict[atom_type][dist[1][k]]).GetIsAromatic()
                                        # print(isAromatic)
                                        if isAromatic == True:
                                            contacts.append(1 / dist[0][k])
                                    '''                            

                dummy_universe.atoms[i].position = (dummy_atom.x, dummy_atom.y, dummy_atom.z)
                    
                dummies_inside.append(i)
                if calculate_bfactor == True:
                    dummy_universe.atoms[i].tempfactor = sum(contacts)
                    #if dummy_atom.is_window == 1:
                    #    dummy_universe.atoms[i].tempfactor = 1.0

    #print("Saving PDB file with the dummy cavity atoms")
    #rdkit.MolToPDBFile(rdkit_molRW, cagePDBout1)
    '''
    if len(dummies_inside)>0:
        for i in range(number_of_dummies):
            if i not in dummies_inside:
                # Lets put all the dummies on the position of the first dummy inside:
                dummy_universe.atoms[i].positions = -np.array([-center_of_mass])
                #dummy_universe.atoms.positions = number_of_dummies * [-center_of_mass]
    else:
        #if therea re not dummies lets put them opposit of orgin
        dummy_universe.atoms.positions = number_of_dummies * [-center_of_mass]
    '''

    if output is not None:
        if output=="tmp/frame":
            MDAnalysis.Merge(dummy_universe.atoms, syst.atoms).atoms.write("{:s}{:05d}".format(output, frame_index))
        else:
            MDAnalysis.Merge(dummy_universe.atoms[dummies_inside], syst.atoms).atoms.write(output)

    # print("Saving MOL file with the dummy cavity atoms")
    # rdkit.MolToMolFile(rdkit_molRW, cageMOLout1)

    #print("--- Total time %s seconds ---" % (time.time() - start_time))

    if len(dummies_inside) > 0:
        cageCavityVolume = cavity_volume(dummy_universe.atoms[dummies_inside].positions, radius=vdwRdummy, volume_grid_size=0.2)
        '''      
        rdkit_molRW = rdkit.RWMol()
        rdkit_conf = rdkit.Conformer(1)
        for dummy_atom in dummy_universe.atoms[dummies_inside].positions:
            indexNewAtom = rdkit_molRW.AddAtom(rdkit.Atom("H"))
            rdkit_conf.SetAtomPosition(indexNewAtom,Point3D(float(dummy_atom[0]),float(dummy_atom[1]),float(dummy_atom[2])))
        rdkit_molRW.AddConformer(rdkit_conf, assignId=True)
        cageCavityVolume = rdkit.ComputeMolVolume(rdkit_molRW, confId=-1, gridSpacing=0.2, boxMargin=2.0)
        '''
        #print("Cage cavity volume = ", cageCavityVolume, " A3")
        
        cavity_hidrophobicity_list = []
        for i, dummy_atom in enumerate(dummy_universe.atoms[dummies_inside]):
            cavity_hidrophobicity_list.append(dummy_universe.atoms[i].tempfactor)
        print("Individual hydrophobicity",cavity_hidrophobicity_list)
        print("Total hydrophobicity", sum(cavity_hidrophobicity_list)/len(cavity_hidrophobicity_list)*cageCavityVolume)
        
        return(cageCavityVolume)
    else:
        print("Cavity with no volume")
        return(0.0)

    # Create an empty editable molecule
    '''
    rdkit_molRW = rdkit.RWMol()
    rdkit_conf = rdkit.Conformer(1)
    for i, dummy_atom in enumerate(calculatedGird.grid):
        if (dummy_atom.vector_angle > 90 and dummy_atom.overlapping_with_cage == 1):
            # if dummy_atom.overlapping_with_cage == 1:
            indexNewAtom = rdkit_molRW.AddAtom(rdkit.Atom(1))
            rdkit_conf.SetAtomPosition(indexNewAtom, Point3D(dummy_atom.x, dummy_atom.y, dummy_atom.z))
            molecInfo = rdkit.AtomPDBResidueInfo()
            molecInfo.SetName(' H')  # the rdkit PDB name has incorrect whitespace
            molecInfo.SetResidueName('HOH')
            molecInfo.SetResidueName(
                ''.ljust(2 - len('H')) + 'HOH')  # the rdkit PDB residue name has incorrect whitespace
            molecInfo.SetResidueNumber(indexNewAtom + 1)
            molecInfo.SetIsHeteroAtom(False)
            molecInfo.SetOccupancy(1.0)
            molecInfo.SetTempFactor(dummy_atom.vector_angle)
            rdkit_molRW.GetAtomWithIdx(indexNewAtom).SetMonomerInfo(molecInfo)

    rdkit_molRW.AddConformer(rdkit_conf, assignId=True)
    print("Saving PDB file with the dummy cavity atoms")
    cagePDBout2 = cagePDBout1.replace(".pdb", "-overlap.pdb")
    rdkit.MolToPDBFile(rdkit_molRW, cagePDBout2)
    '''

    ######################################################################
    ###   Calculate cavity contacts and save to PDB b-factor
    ######################################################################

    # pymol launching: quiet (-q)
    '''
    import pymol
    pymol.pymol_argv = ['pymol', '-q']
    pymol.finish_launching()
    cmd = pymol.cmd
    cmd.load(cageMOL)
    cmd.set('valence', 0)
    cmd.load(cagePDBout1)
    cmd.show("surface", selection=cagePDBout1.replace(".pdb", ""))
    cmd.clip("atoms", 5, "All")
    # def spectrum(expression="count", palette="rainbow",selection="(all)", minimum=None, maximum=None, byres=0, quiet=1):
    cmd.spectrum("b", selection=cagePDBout1.replace(".pdb", ""))
    cmd.load(cagePDBout2)
    cmd.show("spheres", selection=cagePDBout2.replace(".pdb", ""))
    cmd.spectrum("b", selection=cagePDBout2.replace(".pdb", ""))
    cmd.save(cagePDBout1.replace(".pdb", ".pse"))
    '''




import argparse
def get_args():
    parser = argparse.ArgumentParser()
    # Inputs:
    parser.add_argument("-f", default='confout.gro', help="Coordination file (*.pdb, *.gro, *.mol2, *xtc etc.)")
    parser.add_argument("-s", default=None, help="topology file/coordinate file (*gro, *pdb, *tpr etc)")

    # Outputs:
    parser.add_argument("-o", default=None, help="Output coordination file, only pdb supported (due to beta factor)")
    parser.add_argument("-oc", default=None, help="Volumes")

    # Trajectory control:
    parser.add_argument("-n_threads", default=4, help="number of threads for the trajectory calculations")
    parser.add_argument("-b", default=0, help="starting frame")
    parser.add_argument("-e", default=-1, help="last frame")
    parser.add_argument("-stride", default=1, help="stride to calculate trajectory")

    # Parameter control
    parser.add_argument("-grid_spacing", default=1, help="")
    parser.add_argument("-dummy_atom_diameter", default=1.2, help="dummyAtomDiameter")
    parser.add_argument("-distance_threshold_for_90_deg_angle", default=7, help="")
    parser.add_argument("-pore_radius_limit", default=10, help="pore radius limit to consider in cavity")
    parser.add_argument("-save_only_surface", default=True, help="")
    parser.add_argument("-calculate_bfactor", default=True, help="")
    parser.add_argument("-compute_aromatic_contacts", default=False, help="")
    parser.add_argument("-compute_atom_contacts", default=True, help="")
    parser.add_argument("-distThreshold_atom_contacts", default=5.0, help="")
    parser.add_argument("-number_of_dummies", default=500, help="is used for trajectory which requirers constant number of atoms, excess of dummies is put in the center of the cage")
    parser.add_argument("-pymol", default="False", help="open pymol afterwards")
    parser.add_argument("-bcolor", default="False", help="b-factor surface color in pymol")
    parser.add_argument("-bcolor_min", default="None", help="b-factor min value for surface color scale in pymol")
    parser.add_argument("-bcolor_max", default="None", help="b-factor man value for surface color scale in pymol")


    return parser.parse_args()


def single(a,syst, output, grid_spacing = 1.0, distance_threshold_for_90_deg_angle = 7, save_only_surface = True,
         calculate_bfactor = True, compute_aromatic_contacts = False, compute_atom_contacts = True,
         distThreshold_atom_contacts = 5.0, number_of_dummies=500):
    return(a)


if __name__ == '__main__':
    args = get_args()
    if args.s is not None:
        print("Creating temprary folder")
        if args.o is not False:
            os.system("mkdir tmp")
            output= "tmp/frame"
        else:
            output=None
        syst = MDAnalysis.Universe(args.s, args.f)

        run_per_frame = partial(cavity,
                                syst=syst.atoms,
                                output=output,
                                grid_spacing=args.grid_spacing,
                                distance_threshold_for_90_deg_angle=args.distance_threshold_for_90_deg_angle,
                                save_only_surface=args.save_only_surface,
                                calculate_bfactor=args.calculate_bfactor,
                                compute_aromatic_contacts=args.compute_aromatic_contacts,
                                compute_atom_contacts=args.compute_atom_contacts,
                                distThreshold_atom_contacts=args.distThreshold_atom_contacts,
                                number_of_dummies=args.number_of_dummies,
                                pore_radius_limit=args.pore_radius_limit)

        if args.e == -1:
            args.e = syst.trajectory.n_frames


        frame_values = np.arange(int(args.b), int(args.e), int(args.stride))
        #print(frame_values, args.n_threads)


        with Pool(processes=int(args.n_threads)) as pool:
            result = pool.map(run_per_frame, frame_values)
        print(result)

        File = open(args.oc, "w")
        for a in result:
            File.write("{:.2f}\n".format(a))
        File.close()


        FileTraj = MDAnalysis.Writer(args.o)
        for a, ts in enumerate(syst.trajectory):
            syst2 = MDAnalysis.Universe("tmp/frame{:05d}.pdb".format(a))
            temp = MDAnalysis.Merge(syst.atoms, syst2.atoms)
            #temp.write("tmp/file.{:05d}.pdb".format(a))
            FileTraj.write(temp.atoms)
        FileTraj.close()

        os.system("rm -r tmp")

        #pd.DataFrame(results).to_csv("cavity.csv")
    else:
        if args.o is None:
            filePDBout = args.f.replace(".pdb", "_cavity.pdb")
        else:
            filePDBout = args.o
            
        if args.oc is None:
            fileTXTout = args.f.replace(".pdb", "_cavity.txt")
        else:
            fileTXTout = args.oc
        
        cageCavityVolume = cavity(0, MDAnalysis.Universe(args.f), intput = args.f, output = filePDBout,
                                grid_spacing=float(args.grid_spacing),
                                distance_threshold_for_90_deg_angle=float(args.distance_threshold_for_90_deg_angle),
                                distThreshold_atom_contacts=float(args.distThreshold_atom_contacts),
                                number_of_dummies=int(args.number_of_dummies),
                                dummy_atom_diameter=float(args.dummy_atom_diameter),
                                pore_radius_limit=float(args.pore_radius_limit))
        print("Cage cavity volume = ", cageCavityVolume, " A3")
        
        file1 = open(fileTXTout, "w+")
        file1.write("Cage cavity volume = " + "{:.2f}".format(cageCavityVolume) + " A3\n")
        file1.close()

    
    if args.pymol.lower() == "true":
        import pymol
        pymol.pymol_argv = ['pymol', '-q']
        pymol.finish_launching()
        cmd = pymol.cmd
        cmd.load(filePDBout, "cage")
        cmd.set('valence', 0)
        #cmd.set_name(args.o.replace(".pdb", ""), "cage")
        
        cmd.alter('name d', 'vdw="' + str(args.dummy_atom_diameter) + '"')
        
        cmd.create("cavity", selection="elem D", extract=True)
        cmd.show_as("surface", selection="cavity")
        if args.bcolor.lower() == "true":
            if (args.bcolor_min.lower() != "none" and args.bcolor_max.lower() != "none"):
                cmd.spectrum("b", minimum=float(args.bcolor_min), maximum=float(args.bcolor_max), selection="cavity")
                cmd.ramp_new("ramp", "cavity", [-float(args.bcolor_min),float(args.bcolor_max)], "rainbow")
                # def spectrum(expression="count", palette="rainbow",selection="(all)", minimum=None, maximum=None, byres=0, quiet=1):

            else:
                cmd.spectrum("b", selection="cavity")

        else:
            cmd.set("surface_color", "cyan", selection="cavity")
            cmd.set("transparency", 0.5, selection="cavity")
        
        cmd.show("surface", selection="cage")
        cmd.set("surface_color", "green", selection="cage") 
        cmd.set("transparency", 0.7, selection="cage")
        cmd.hide("surface", selection="cage")
        
        cmd.clip("atoms", 5, "All")
        
        cmd.orient("cage")
        cmd.zoom("cage")
   
        
        #cmd.load(cagePDBout2)
        #cmd.show("spheres", selection=cagePDBout2.replace(".pdb", ""))
        #cmd.spectrum("b", selection=cagePDBout2.replace(".pdb", ""))
        cmd.save(filePDBout.replace(".pdb", ".pse"))
