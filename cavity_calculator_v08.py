import os
import rdkit.Chem.AllChem as rdkit
import rdkit.Chem.rdmolops as rdkitrdmolops
from rdkit.Geometry import Point3D
from rdkit.Chem import Descriptors
from copy import deepcopy
import numpy as np
import math
import time
from scipy.spatial import KDTree

########## cavity calculation
cageMOL = "cage_ACIE_2006_45_901.mol2" #we need mol files as pdb file was not behaving OK with the aromatic
cagePDB = cageMOL.replace(".mol2", ".pdb")
print(cagePDB)
cagePDBout1 = cagePDB.replace(".pdb", ".cavity.pdb")

grid_spacing = 1.0
distance_threshold_for_90_deg_angle = 7 #10
save_only_surface = True
calculate_bfactor = True
compute_aromatic_contacts = False
compute_atom_contacts = True
distThreshold_atom_contacts = 5.0

calculate_windows = True

threads_KDThree = 4

######################################################################
###   Calculate cavity
######################################################################
print("\n--- Calculation of the cavity ---")

start_time = time.time()

## Using the PDB file as the MOL file for this example is not working, some issues as it was generated from the CIF file
#rdkit_cage = rdkit.MolFromMol2File(cageMOL,removeHs=False)
rdkit_cage = rdkit.MolFromPDBFile(cagePDB,sanitize=False,removeHs=False,proximityBonding=False)


#Calculate the centeor of mass
atoms = [ atom for atom in rdkit_cage.GetAtoms() ]
numberOfAtoms = rdkit_cage.GetNumAtoms()
xyzPositionArray =  np.array([list( rdkit_cage.GetConformer().GetAtomPosition(i) ) for i in range(numberOfAtoms)])
center_of_mass = np.array(sum(atoms[i].GetMass()*xyzPositionArray[i] for i in range(numberOfAtoms)))/Descriptors.MolWt(rdkit_cage)
print("center_of_mass= ", center_of_mass)


# Get the xyx atom coordinates
xyzAtoms = np.array([rdkit_cage.GetConformer().GetAtomPosition(cage_atom.GetIdx()) for cage_atom in rdkit_cage.GetAtoms()], dtype=np.float64)
# Calculate neighbors with KDTree. Leaf_size affects the speed of a query and the memory required to store the constructed tree. 
# The amount of memory needed to store the tree scales as approximately n_samples / leaf_size.
kdtxyzAtoms = KDTree(xyzAtoms, leafsize = 20)
#Calculate the set of atoms within the euclidean distanc upper bound
result = kdtxyzAtoms.query(xyzAtoms[30], k = None, p = 2, distance_upper_bound = 3.0)

#We use the center of mass of the molecule as the pore center of mass and the pore radius as the 95% of the close contact to the center_of_mass
pore_center_of_mass = center_of_mass
distancesFromCOM = kdtxyzAtoms.query(center_of_mass, k = None, p = 2)
pore_radius = distancesFromCOM[0][0]*0.8
maxDistanceFromCOM = (distancesFromCOM[0][-1])
pore_diameter = pore_radius*2

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
        
        self.d_from_center = np.linalg.norm( np.array(self.pos) -  np.array(center) )
        
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
        self.n = 2*self.points + 1  
        self.grid = [GridPoint(i,j,k, self.points, self.center, self.grid_spacing) for k in range(self.n) for j in range(self.n) for i in range(self.n)]
        self.grid = np.array(self.grid)
        self.gridPosList = []
        for i in self.grid:
            self.gridPosList.append(i.pos)


box_size = maxDistanceFromCOM
print("box_size=", box_size)
calculatedGird = CageGrid(pore_center_of_mass, box_size, delta = 0, grid_spacing = grid_spacing)

#-#-#-#-# KDTree algorithm to speed up the method, needs to be implemented with the code below

# Create a KDThree of the dummy atom gird
calculatedGirdKDTree = KDTree(calculatedGird.gridPosList, leafsize = 20)
# Calculate contacts of dummy atoms with a distance < 1.1*grid_spacing, rages from 1 to 7 (as 1 contact is always, the atom itself)

'''
for dummy_atom in calculatedGird.grid:
    xyzDummySet = calculatedGirdKDTree.query(dummy_atom.pos, k=None, p=2, distance_upper_bound=1.1*grid_spacing)
    dummy_atom.neighbors = xyzDummySet[1] #does not work as the index does not correspond
'''

## IMPROVED ALGORITHM USING KDThree
atom_type_list = []
coords_dict = {}
vdwR_dict = {}
atom_idx_dict = {}

for cage_atom in rdkit_cage.GetAtoms():
    atom_type = cage_atom.GetSymbol()
    atom_idx = cage_atom.GetIdx()
    pos = rdkit_cage.GetConformer().GetAtomPosition(atom_idx)

    if atom_type not in atom_type_list:
        atom_type_list.append(atom_type)
        vdwR = rdkit.GetPeriodicTable().GetRvdw(cage_atom.GetAtomicNum())
        vdwR_dict.setdefault(atom_type, []).append(vdwR)
        
    coords_dict.setdefault(atom_type, []).append(list(pos))
    atom_idx_dict.setdefault(atom_type, []).append(atom_idx)

#print(atom_type_list)
print(vdwR_dict)
#print(coords_dict)
#print(atom_idx_dict)

#Create a KDTree for each atom type and store in a dict
KDTree_dict = {}
for atom_type in atom_type_list:
    #print(coords_dict[atom_type])
    #print(atom_type)
    KDTree_dict[atom_type] = KDTree(coords_dict[atom_type], leafsize = 20)

#Radius of the dummy H atoms in the cavity
vdwRdummy = rdkit.GetPeriodicTable().GetRvdw(1)
for i in calculatedGird.grid:
    dist1 = np.linalg.norm(np.array(pore_center_of_mass)-np.array(i.pos))
    if dist1 > 0:
            vect1Norm = (np.array(pore_center_of_mass)-np.array(i.pos)) / dist1
    if dist1 < pore_radius:        
            i.inside_cavity = 1      
    if i.inside_cavity == 0:
        summAngles = []
        distancesAngles = []      
        for atom_type in atom_type_list:
            xyzAtomsSet2 = KDTree_dict[atom_type].query(i.pos, k = None, p = 2, distance_upper_bound = distance_threshold_for_90_deg_angle, workers = threads_KDThree)
            if xyzAtomsSet2[1]:
                vdwR = vdwR_dict[atom_type][0]
                distThreshold = vdwR + vdwRdummy
                if xyzAtomsSet2[0][0] < distThreshold:
                    #print("Dummy atom overlapping with cage")
                    i.overlapping_with_cage = 1
                    #break
                '''
                for z in range(0, len(xyzAtomsSet2[1])):
                    #print(xyzAtomsSet2[0][0], distThreshold, vdwR, vdwRdummy)
                    if xyzAtomsSet2[0][z] < 1.01*distThreshold:
                        #print("Dummy atom overlapping with cage")
                        i.overlapping_with_cage = 1
                        break
                '''
                        
                        
                # check that distances are the same by KDTree query an position based calculation
                """"
                for z in range (0,len(xyzAtomsSet2[1])):
                    atom_pos = coords_dict[atom_type][xyzAtomsSet2[1][z]]
                    dist2 = np.linalg.norm(np.array(atom_pos)-np.array(i.pos))
                    print(xyzAtomsSet2[0][z], dist2)
                """
                
                for atom_pos_index in xyzAtomsSet2[1]:
                    atom_pos = coords_dict[atom_type][atom_pos_index]
                    dist2 = np.linalg.norm(np.array(atom_pos)-np.array(i.pos))
                    vect2Norm = (np.array(atom_pos)-np.array(i.pos)) / dist2
                    angle = np.arccos(np.dot(vect1Norm, vect2Norm))
                    summAngles.append(angle)
                    distancesAngles.append(1/dist2)
        if (summAngles):
            #averageSummAngles = sum(summAngles) / len(summAngles)
            averageSummAngles = np.average(summAngles, axis=None, weights=distancesAngles)
            averageSummAngles_deg = np.degrees(averageSummAngles)
            i.vector_angle = averageSummAngles_deg 
            if i.overlapping_with_cage == 0:
                if averageSummAngles_deg > 90:
                    #print("Added dummy atom to cavity", np.degrees(averageSummAngles), i.overlapping_with_cage)
                    i.inside_cavity = 1

      
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
for i in calculatedGird.grid:
    if i.inside_cavity == 1:
        calculatedGirdContacts.append(i.pos) 
calculatedGirdContactsKDTree = KDTree(calculatedGirdContacts, leafsize = 20)

# Calculate contacts of dummy atoms with a distance < 1.1*grid_spacing, rages from 1 to 7 (as 1 contact is alwas, the atom itself)
for i, dummy_atom in enumerate(calculatedGird.grid):
    if dummy_atom.inside_cavity == 1:
        xyzDummySet = calculatedGirdContactsKDTree.query(dummy_atom.pos, k=None, p=2, distance_upper_bound=1.1*grid_spacing)
        dummy_atom.number_of_neighbors = len(xyzDummySet[1])


# Create a KDThree of the dummy atom girds to that overlap with the cage

for i in calculatedGird.grid:   
    for atom_type in atom_type_list:
        vdwR = vdwR_dict[atom_type][0]
        distThreshold = vdwR + vdwRdummy
        xyzAtomsSet2 = KDTree_dict[atom_type].query(i.pos, k = None, p = 2, distance_upper_bound = distThreshold, workers = threads_KDThree)
        if xyzAtomsSet2[0]:
            i.overlapping_with_cage = 1
        else:
            i.overlapping_with_cage = 0


overlappingCalculatedGirdContacts = []
overlappingCalculatedGirdContacts_angles = []
for i in calculatedGird.grid:
    if i.overlapping_with_cage == 1:
        #if i.vector_angle > 90:
        overlappingCalculatedGirdContacts.append(i.pos)
        overlappingCalculatedGirdContacts_angles.append(i.vector_angle)
overlappingCalculatedGirdContactsKDTree = KDTree(overlappingCalculatedGirdContacts, leafsize = 20)

## IMPROVED ALGORITHM USING KDThree
# Remove dummy atoms overapping by the cage

#distance111 = KDTree_dict['Pd'].query([0,0,0], k = None, p = 2, distance_upper_bound = 50)
#print(distance111)
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


#Create an empty editable molecule
rdkit_molRW = rdkit.RWMol()
rdkit_conf = rdkit.Conformer(1)

for i, dummy_atom in enumerate(calculatedGird.grid):
    if (dummy_atom.inside_cavity == 1 and dummy_atom.overlapping_with_cage == 0):
        if calculate_windows == True:
            distWindow = overlappingCalculatedGirdContactsKDTree.query(dummy_atom.pos, k = None, p = 2, distance_upper_bound = 1.1*grid_spacing, workers = threads_KDThree)
            #print(distWindow)
            if len(distWindow[0]) == 0:
                dummy_atom.is_window = 1
            if len(distWindow[0]) > 0:
                if overlappingCalculatedGirdContacts_angles[distWindow[1][0]] < 100: ## 100 deg is necessary to include the atoms close to the surface
                    dummy_atom.is_window = 1
                    #print("dummy atom not in a window")

        # Remove isolated dummy atoms with a distance < 1.1*grid_spacing
        if (  (dummy_atom.number_of_neighbors > 1 and dummy_atom.number_of_neighbors < 7 and save_only_surface == True)):
            if calculate_bfactor == True:
                contacts = []                
                for atom_type in atom_type_list:
                    if atom_type != "H":
                        dist = KDTree_dict[atom_type].query(dummy_atom.pos, k = None, p = 2, distance_upper_bound = distThreshold_atom_contacts, workers = threads_KDThree)

                        if dist[0]:
                            for k in range (0, len(dist[0])):
                                if compute_atom_contacts == True:
                                    #print("compute_atom_contacts")
                                    contacts.append(1/dist[0][k])
                                elif compute_aromatic_contacts == True:
                                    isAromatic = rdkit_cage.GetAtomWithIdx(atom_idx_dict[atom_type][dist[1][k]]).GetIsAromatic()
                                    #print(isAromatic)
                                    if isAromatic == True:
                                        contacts.append(1/dist[0][k])
                                            
            indexNewAtom = rdkit_molRW.AddAtom(rdkit.Atom(1))
            #Change the postion of the atom to a custom x, y, z postion
            rdkit_conf.SetAtomPosition(indexNewAtom,Point3D(dummy_atom.x,dummy_atom.y,dummy_atom.z))
            molecInfo = rdkit.AtomPDBResidueInfo()            
            molecInfo.SetName(' H') # the rdkit PDB name has incorrect whitespace
            molecInfo.SetResidueName('HOH')            
            molecInfo.SetResidueName(''.ljust(2-len('H'))+'HOH') # the rdkit PDB residue name has incorrect whitespace
            molecInfo.SetResidueNumber(indexNewAtom+1)
            molecInfo.SetIsHeteroAtom(False)
            molecInfo.SetOccupancy(1.0)
            if calculate_bfactor == True:
                #molecInfo.SetTempFactor( sum(contacts) )  #Save in the b-factor the mumber of contacts
                
                if dummy_atom.is_window == 1:
                    molecInfo.SetTempFactor( 0 )  #Save in the b-factor the mumber of contacts
                else:
                    molecInfo.SetTempFactor( 1 )
            else:
                molecInfo.SetTempFactor( 1 )
            #molecInfo.SetTempFactor(dummy_atom.inside_cavity)
            #molecInfo.SetTempFactor(dummy_atom.d_from_center) # Save in the b-factor the distance from the center
            rdkit_molRW.GetAtomWithIdx(indexNewAtom).SetMonomerInfo(molecInfo)

rdkit_molRW.AddConformer(rdkit_conf, assignId=True)


print("Saving PDB file with the dummy cavity atoms")
rdkit.MolToPDBFile(rdkit_molRW, cagePDBout1)
#print("Saving MOL file with the dummy cavity atoms")
#rdkit.MolToMolFile(rdkit_molRW, cageMOLout1)

if len(rdkit_molRW.GetAtoms()) >0:
	cageCavityVolume = rdkit.ComputeMolVolume(rdkit_molRW, confId=-1, gridSpacing=0.2, boxMargin=2.0)
	print("Cage cavity volume = ", cageCavityVolume, " A3")
	cage_name = cagePDB.replace(".pdb", "")
	file1 = open(cage_name+"_cavity_vol_calc.txt","w+")
	file1.write("PCage cavity volume = " + str(cageCavityVolume) + " A3\n")
	file1.close()
else:
	print("Cavity with no volume")





#Create an empty editable molecule
rdkit_molRW = rdkit.RWMol()
rdkit_conf = rdkit.Conformer(1)
for i, dummy_atom in enumerate(calculatedGird.grid):
    if (dummy_atom.vector_angle > 90 and dummy_atom.overlapping_with_cage == 1):
    #if dummy_atom.overlapping_with_cage == 1:                    
            indexNewAtom = rdkit_molRW.AddAtom(rdkit.Atom(1))
            rdkit_conf.SetAtomPosition(indexNewAtom,Point3D(dummy_atom.x,dummy_atom.y,dummy_atom.z))
            molecInfo = rdkit.AtomPDBResidueInfo()            
            molecInfo.SetName(' H') # the rdkit PDB name has incorrect whitespace
            molecInfo.SetResidueName('HOH')            
            molecInfo.SetResidueName(''.ljust(2-len('H'))+'HOH') # the rdkit PDB residue name has incorrect whitespace
            molecInfo.SetResidueNumber(indexNewAtom+1)
            molecInfo.SetIsHeteroAtom(False)
            molecInfo.SetOccupancy(1.0)
            molecInfo.SetTempFactor( dummy_atom.vector_angle )
            rdkit_molRW.GetAtomWithIdx(indexNewAtom).SetMonomerInfo(molecInfo)

rdkit_molRW.AddConformer(rdkit_conf, assignId=True)
print("Saving PDB file with the dummy cavity atoms")
cagePDBout2 = cagePDBout1.replace(".pdb", "-overlap.pdb")
rdkit.MolToPDBFile(rdkit_molRW, cagePDBout2)






print("--- Total time %s seconds ---" % (time.time() - start_time))

######################################################################
###   Calculate cavity contacts and save to PDB b-factor
######################################################################

# pymol launching: quiet (-q)
import pymol
pymol.pymol_argv = ['pymol','-q']
pymol.finish_launching()
cmd = pymol.cmd
#cmd.load(cageMOL)
cmd.load(cagePDB)
cmd.set('valence', 0)
cmd.load(cagePDBout1)
cmd.show("surface", selection=cagePDBout1.replace(".pdb", ""))
cmd.clip("atoms", 5, "All")
#def spectrum(expression="count", palette="rainbow",selection="(all)", minimum=None, maximum=None, byres=0, quiet=1):
cmd.spectrum("b", selection=cagePDBout1.replace(".pdb", ""))
cmd.load(cagePDBout2)
cmd.show("spheres", selection=cagePDBout2.replace(".pdb", ""))
cmd.spectrum("b", selection=cagePDBout2.replace(".pdb", ""))

cage_name = cagePDB.replace(".pdb", "")
cmd.show("surface", cage_name)
cmd.set("surface_color", "withe", cage_name)
#cmd.set("transparency", 0.5, cage_name)

cmd.save(cagePDBout1.replace(".pdb", ".pse"))
