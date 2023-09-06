import numpy as np
import re

from tempfile import mkdtemp
from rdkit import Chem
import rdkit.Chem.AllChem as rdkit
from rdkit.Geometry import Point3D

from CageCavityCalc.data import atom_mass, vdw_radii, name_to_atomic_number
from CageCavityCalc.log import logger

def atom_names_to_masses(names):
    atom_masses = []
    for name in names:
        if name in atom_mass:
            atom_masses.append(atom_mass[name])
        else:
            atom_masses.append(0.0)
    return np.array(atom_masses)


def atom_names_to_vdw(names):
    atom_vdw_radii = []
    for name in names:
        if name in vdw_radii:
            atom_vdw_radii.append(vdw_radii[name])
        else:
            atom_vdw_radii.append(0.0)
    return np.array(atom_vdw_radii)


# ----------------- INPUT -------------------------


def read_positions_and_atom_names_from_file(filename):
    positions = None
    atom_names = None
    if filename.endswith(".pdb"):
        positions, atom_names = read_pdb(filename)
    elif filename.endswith(".mol2"):
        positions, atom_names = read_mol2(filename)
    else:
        positions, atom_names = read_other(filename)

    return positions, atom_names, atom_names_to_masses(atom_names), atom_names_to_vdw(atom_names)


def read_pdb(filename):
    positions = []
    atom_names = []
    with open(filename) as File:
        text = File.read()
        for line in text.splitlines():
            if line.split()[0] == "HETATM" or line.split()[0] == "ATOM":
                temp = np.array(list(map(float, line[30:54].split())))
                positions.append(temp)
                name_and_number = line[12:16].upper()
                name_strip_number = re.match('\s*([A-Z]+)', name_and_number).group(1)
                atom_names.append(name_strip_number)
            elif line.split()[0] == "END":
                break

    return np.array(positions), atom_names


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
                name_and_number = line.split()[1].upper()
                name_strip_number = re.match('([A-Z]+)', name_and_number).group(1)
                atom_names.append(name_strip_number)

    return np.array(positions), atom_names


def read_other(filename):
    try:
        import MDAnalysis
    except:
        print("The other formats are supported by MDAnalysis, which has been not found")
        exit()
    syst = MDAnalysis.Universe(filename)

    # we take only first two inputs (positions and atoms)
    return read_mdanalysis(syst)[:2]


def read_cgbind(cgbind_cage):
    try:
        import cgbind
    except:
        print("Could not load cgbind")
        exit()

    atom_names = [atom.label.upper() for atom in cgbind_cage.atoms]
    positions = cgbind_cage.get_coords()
    return np.array(positions), atom_names, atom_names_to_masses(atom_names), atom_names_to_vdw(atom_names)


def read_mdanalysis(syst):
    try:
        import MDAnalysis
    except:
        print("Could not load MDAnalysis")
        exit()

    atom_names = []
    for name in syst.atoms.names:
        name_strip_number = re.match('([A-Z]+)', name.upper()).group(1)
        atom_names.append(name_strip_number)
    positions = syst.atoms.positions

    return np.array(positions), atom_names, atom_names_to_masses(atom_names), atom_names_to_vdw(atom_names)

def read_positions_and_atom_names_from_array(positions, atom_names):

    just_atom_names = [re.match('([A-Z]+)', name_and_number.upper()).group(1) for name_and_number in atom_names]
    return positions, just_atom_names, atom_names_to_masses(just_atom_names), atom_names_to_vdw(just_atom_names)


# ----------------- OUTPUT -------------------------


def print_to_file(filename, positions, atom_names, property_values = None):

    if property_values is not None and not filename.endswith('.pdb'):
        print("properties can be saved only to pdb file!")
        exit()

    if filename.endswith(".xyz"):
        print_to_xyz_file(filename, positions, atom_names)
    elif filename.endswith(".pdb"):
        print_to_pdb_file(filename, positions, atom_names, property_values)
    else:
        print_to_other_file(filename, positions, atom_names)  # TODO


def print_to_pdb_file(filename, positions, atom_names, property_values=None):
    with open(filename, 'w') as xyz_file:
        # print(len(positions), "CageCavityCalc", sep='\n', file=xyz_file)
        for a, pos in enumerate(positions):

            element = re.search("([a-zA-Z]+)", atom_names[a]).group(1)

            if atom_names[a] != "D":
                print(
                    f"ATOM  {a:>5d} {atom_names[a].upper():<4s}  CG A   0    {pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}{0:6.2f}{0.0:6.2f}          {element.upper():>2s}",
                    file=xyz_file)
            else:
                if property_values is not None:
                    print(f"ATOM  {a:>5d} {atom_names[a].upper():<4s}  CV B   1    {pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}{0:6.2f}{property_values[a]:6.2f}          {element.upper():>2s}",file=xyz_file)
                else:
                    print(
                        f"ATOM  {a:>5d} {atom_names[a].upper():<4s}  CV B   1    {pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}{0:6.2f}{0:6.2f}          {element.upper():>2s}", file=xyz_file)
    return None


def print_to_xyz_file(filename, positions, atom_names):
    """
    Print a standard .xyz file from a set of atoms
    :param filename: (str)
    :param positions:
    :param atom_names:
    """

    with open(filename, 'w') as xyz_file:
        print(len(positions), "CageCavityCalc", sep='\n', file=xyz_file)
        for a, pos in enumerate(positions):
            print(f'{atom_names[a].upper():<3} {pos[0]:^10.5f} {pos[1]:^10.5f} {pos[2]:^10.5f}', file=xyz_file)
    return None


def print_to_other_file(filename, positions, atom_names):
    try:
        import MDAnalysis
    except:
        print("The other formats than *.pdb and *.xyz are supported by MDAnalysis, which has been not found")
        exit()

    n_atoms = len(positions)


    resids = [1 * (atom_name == 'D') for atom_name in atom_names]
    resnames = ['CG'*(atom_name!='D')+'CV'*(atom_name=='D') for atom_name in atom_names]

    n_resids=len(set(resids))

    new_universe = MDAnalysis.Universe.empty(n_atoms, n_residues=n_resids, atom_resindex=resids, trajectory=True)
    new_universe.add_TopologyAttr('tempfactors')

    # we put dummies far away from cage:
    new_universe.atoms.positions = positions
    new_universe.add_TopologyAttr('name', atom_names)

    new_universe.add_TopologyAttr('resname', ['CG',  'CV'][:n_resids])
    new_universe.atoms.write(filename)
    return None


def convert_to_mol2(positions, atom_names):
    # openbabel
    try:
        from openbabel import openbabel
    except:
        print("You must either provide .mol2 file or install openbabel")


    tmpdir_path = mkdtemp()
    print_to_pdb_file(tmpdir_path+"/temp.pdb", positions, atom_names)

    # Convert using Open Babel to .mol2 file
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("pdb", "mol2")
    mol = openbabel.OBMol()
    ob_conversion.ReadFile(mol, tmpdir_path + "/temp.pdb")
    ob_conversion.WriteFile(mol, tmpdir_path + "/temp.mol2")

    rdkit_cage = rdkit.MolFromMol2File(tmpdir_path + "/temp.mol2", removeHs=False)
    return rdkit_cage

def print_pymol_file(filename, property_values=None, dummy_atom_radii=1):
    if not filename.endswith(".pml"):
        logger.warning("To easy open this file in pymol it should have *.pml extension")


    with open(filename, 'w') as file:
        print(f"load {filename[:filename.find('.')]}.pdb", file=file)
        print("extract cavity, resname CV", file=file)
        print("alter name D, vdw="+str(dummy_atom_radii), file=file)
        print("show_as surface, cavity", file=file)

        if property_values is not None:
            print(f"spectrum b, blue_white_red,cavity, minimum={np.min(property_values):f}, maximum={np.max(property_values):f}", file=file)
            print(f"ramp_new 'ramp', cavity, [{np.min(property_values):f},{np.mean(property_values):f},{np.max(property_values):f}], ['blue','white','red']", file=file)
            print("recolor", file=file)


        ######################################################################
        ###        Open the saved cavity in PyMol
        ######################################################################


def run_pymol():
        # pymol launching: quiet (-q)
        import pymol
        pymol.pymol_argv = ['pymol','-q']
        pymol.finish_launching()
        cmd = pymol.cmd
        cmd.load(cageMOL, "cage")
        cmd.set('valence', 0)
        cmd.load(cagePDBout1, "cavity")
        if(printLevel == 2):
            cmd.load(cagePDB.replace(".pdb", "_box_cavity.pdb"), "box")

        cmd.alter('name D', 'vdw="' + str(dummy_atom_radii) + '"')

        #cmd.show_as("nonbonded", selection="cavity")
        #cmd.show_as("surface", selection="cavity")
        cmd.show_as("spheres", selection="cavity")
        #cmd.show_as("spheres", selection="cage")

        cmd.spectrum("b", selection="cavity",palette="blue_white_red",minimum=min(hydrophobicity_cavity_dummy_atoms), maximum=max(hydrophobicity_cavity_dummy_atoms))
        #cmd.set("surface_color", "withe", selection="cavity")
        #cmd.set("transparency", 0.5, cage_name)

        cmd.ramp_new("ramp", "cavity", [min(hydrophobicity_cavity_dummy_atoms),(min(hydrophobicity_cavity_dummy_atoms)+max(hydrophobicity_cavity_dummy_atoms))/2,max(hydrophobicity_cavity_dummy_atoms)], ["blue","white","red"] )
        cmd.recolor()

        cmd.clip("atoms", 5, "All")
        cmd.orient("cage")
        cmd.zoom("cage")

        cmd.save(cagePDBout1.replace(".pdb", ".pse"))
