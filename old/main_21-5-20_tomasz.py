import numpy as np
import math
import time
from scipy.spatial import KDTree
import MDAnalysis
import os
from functools import partial
from multiprocess import Pool

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

def cavity(frame_index, syst, output, grid_spacing = 1.0, distance_threshold_for_90_deg_angle = 7, save_only_surface = True,
         calculate_bfactor = True, compute_aromatic_contacts = False, compute_atom_contacts = True,
         distThreshold_atom_contacts = 5.0, number_of_dummies=500):

    syst.universe.trajectory[frame_index]
    if frame_index>0:
        print("Progress: {:d}/{:d}".format(frame_index, len(syst.universe.trajectory)), end = '\r')

    ########## cavity calculation
    cageMOL = "fdsfdsfdsf.pdb"
    # we need mol files as pdb file was not behaving OK with the aromatic
    cagePDB = cageMOL.replace(".mol2", ".pdb")
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
    #rdkit_cage = rdkit.MolFromMol2File(cageMOL, removeHs=False)
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
    vdwRdummy = vdw_radii['h']# rdkit.GetPeriodicTable().GetRvdw(1)
    for a, i in enumerate(calculatedGird.grid):
        dist1 = np.linalg.norm(np.array(pore_center_of_mass) - np.array(i.pos))
        if dist1 > 0:
            vect1Norm = (np.array(pore_center_of_mass) - np.array(i.pos)) / dist1
        if dist1 < pore_radius:
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
            if (
            (dummy_atom.number_of_neighbors > 1 and dummy_atom.number_of_neighbors < 7 and save_only_surface == True)):
                if calculate_bfactor == True:
                    contacts = []
                    for atom_type in atom_type_list:
                        if atom_type != "H":
                            dist = KDTree_dict[atom_type].query(dummy_atom.pos, k=None, p=2,
                                                                distance_upper_bound=distThreshold_atom_contacts)

                            if dist[0]:
                                for k in range(0, len(dist[0])):
                                    if compute_atom_contacts == True:
                                        # print("compute_atom_contacts")
                                        contacts.append(1 / dist[0][k])
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
                    if dummy_atom.is_window == 1:
                        dummy_universe.atoms[i].tempfactor = 1.0

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
        cageCavityVolume = cavity_volume(dummy_universe.atoms[dummies_inside].positions, vdwRdummy)
        print("Cage cavity volume = ", cageCavityVolume, " A3")
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
    parser.add_argument("-oc", default='cavity.dat', help="Volumes")

    # Trajectory control:
    parser.add_argument("-n_threads", default=4, help="number of threads for the trajectory calculations")
    parser.add_argument("-b", default=0, help="starting frame")
    parser.add_argument("-e", default=-1, help="last frame")
    parser.add_argument("-stride", default=1, help="stride to calculate trajectory")

    # Parameter control
    parser.add_argument("-grid_spacing", default=1, help="")
    parser.add_argument("-distance_threshold_for_90_deg_angle", default=7, help="")
    parser.add_argument("-save_only_surface", default=True, help="")
    parser.add_argument("-calculate_bfactor", default=True, help="")
    parser.add_argument("-compute_aromatic_contacts", default=False, help="")
    parser.add_argument("-compute_atom_contacts", default=True, help="")
    parser.add_argument("-distThreshold_atom_contacts", default=5.0, help="")
    parser.add_argument("-number_of_dummies", default=500, help="is used for trajectory which requirers constant number of atoms, excess of dummies is put in the center of the cage")
    parser.add_argument("-pymol", default=False, help="open pymol afterwards")


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
                                number_of_dummies=args.number_of_dummies)

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
        cageCavityVolume = cavity(0, MDAnalysis.Universe(args.f), output = args.o)
        print("Cage cavity volume = ", cageCavityVolume, " A3")

    if args.pymol==True:
        import pymol
        pymol.pymol_argv = ['pymol', '-q']
        pymol.finish_launching()
        cmd = pymol.cmd
        cmd.load(args.o)
        cmd.set('valence', 0)
        cmd.show("surface", selection="elem D")
        cmd.clip("atoms", 5, "All")
        # def spectrum(expression="count", palette="rainbow",selection="(all)", minimum=None, maximum=None, byres=0, quiet=1):
        cmd.spectrum("b", selection="elem D")
        #cmd.load(cagePDBout2)
        #cmd.show("spheres", selection=cagePDBout2.replace(".pdb", ""))
        #cmd.spectrum("b", selection=cagePDBout2.replace(".pdb", ""))
        #cmd.save(cagePDBout1.replace(".pdb", ".pse"))
