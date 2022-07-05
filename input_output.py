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

