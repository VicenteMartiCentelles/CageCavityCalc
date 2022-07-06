import numpy as np
import re
from data import atom_mass, vdw_radii

def read_pdb():
    None

def read_mol2(filename):
    positions = []
    atom_names = []
    with open(filename) as File:
        text = File.read()

        text = text[text.find("@<TRIPOS>ATOM") + 14:]
        text = text[:text.find("@<TRIPOS>")]

        for line in text.splitlines():
            if len(line) > 0:
                positions.append([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])
                name_and_number = line.split()[1].lower()
                name_strip_number = re.match('([a-z]+)', name_and_number).group(1)
                atom_names.append(name_strip_number)

    return np.array(positions), atom_names
'''
self.positions = None
self.atom_names = None
self.atom_masses = None
self.atom_vdw = None
self.n_atoms = 0
'''
def atom_names_to_masses(names):
    atom_masses = []
    for name in names:
        atom_masses.append(atom_mass[name])
    return np.array(atom_masses)

def atom_names_to_vdw(names):
    atom_vdw_radii = []
    for name in names:
        atom_vdw_radii.append(vdw_radii[name])
    return np.array(atom_vdw_radii)


def read_positions_and_atom_names_from_file(filename):
    positions = None
    atom_names = None
    if filename.endswith(".pdb"):
        positions, atom_names = read_pdb(filename)
    elif filename.endswith(".mol2"):
        positions, atom_names = read_mol2(filename)
    else:
        positions, atom_names = read_other(filename) #TODO

    return positions, atom_names,  atom_names_to_masses(atom_names), atom_names_to_vdw(atom_names)


def print_to_file(filename, positions, atom_names):
    if filename.endswith(".xyz"):
        print_to_xyz_file(filename, positions, atom_names)
    elif filename.endswith(".pdb"):
        print_to_pdb_file(filename, positions, atom_names)


def print_to_pdb_file(filename, positions, atom_names):
    with open(filename, 'w') as xyz_file:
        #print(len(positions), "CageCavityCalc", sep='\n', file=xyz_file)
        for a, pos in enumerate(positions):
            if atom_names[a]!="D":
                print(f"ATOM  {a:>5d} {atom_names[a].upper()+str(a):<4s}  CG A   0    {pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}{0:6.2f}{0:6.2f}", file=xyz_file)
            else:
                print(f"ATOM  {a:>5d} {atom_names[a].upper():<4s}  CV B   1    {pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}{0:6.2f}{0:6.2f}", file=xyz_file)
    return None
'''
print(f"{'ATOM':>5s}{a:>12d}{atom_names[a].upper()+str(a):>7s}{1:>8d}{:>6s}{:>5s}{:>12.6f}{:>11.6f}")
print("ATOM   {a:>5d}  {atom_names[a]:>4s} CG A    0     {pos[0]:>8.3f}{pos[0]:>8.3f}{pos[0]:>8.3f} {0:6.2f}{0:6.2f}"
print("ATOM   {a:>5d}  {atom_names[a]:>4s} CV B    1     {pos[0]:>8.3f}{pos[0]:>8.3f}{pos[0]:>8.3f} {0:6.2f}{0:0.2f}"  

ATOM     24  C1  DMSO   78      20.530  14.660  20.390  1.00  0.00
ATOM    134  H134  CG A  0      -7.407  -7.407  -7.407  0.00  0.00
ATOM    309  DM2 UNL M   1      10.457  11.039  31.690  0.00  1.00      M    D

ATOM 	1-4	“ATOM”		character
7-11#	Atom serial number	right	integer
13-16	Atom name	left*	character
17	Alternate location indicator		character
18-20§	Residue name	right	character
22	Chain identifier		character
23-26	Residue sequence number	right	integer
27	Code for insertions of residues		character
31-38	X orthogonal Å coordinate	right	real (8.3)
39-46	Y orthogonal Å coordinate	right	real (8.3)
47-54	Z orthogonal Å coordinate	right	real (8.3)
55-60	Occupancy	right	real (6.2)
61-66	Temperature factor	right	real (6.2)
73-76	Segment identifier¶	left	character
77-78	Element symbol	right	character

print("{:>5s}{:>12s}{:>7s}{:>8s}{:>6s}{:>5s}{:>12.6f}{:>11.6f}".format(
    modified_line[0],
    modified_line[1],
    modified_line[2],
    name,
    modified_line[4],
    modified_line[5],
    float(modified_line[6]),
    float(modified_line[7]),
),
'''
def print_to_xyz_file(filename, positions, atom_names):
    """
    Print a standard .xyz file from a set of atoms
    :param filename: (str)
    :param positions:
    :param atom_names:
    """
    print(positions.shape, atom_names.shape)
    print(positions)
    print(atom_names)
    with open(filename, 'w') as xyz_file:
        print(len(positions), "CageCavityCalc", sep='\n', file=xyz_file)
        for a, pos in enumerate(positions):
            print(f'{atom_names[a].upper():<3} {pos[0]:^10.5f} {pos[1]:^10.5f} {pos[2]:^10.5f}', file=xyz_file)
    return None



