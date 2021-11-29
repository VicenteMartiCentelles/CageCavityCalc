import os
import rdkit.Chem.AllChem as rdkit
import rdkit.Chem.rdmolops as rdkitrdmolops
from rdkit.Geometry import Point3D
from rdkit.Chem import Descriptors
from rdkit import Chem
from copy import deepcopy
import numpy as np
import math
import time
from scipy.spatial import KDTree

from hydrophobicity_values import hydrophValues
from calculations import cavity_volume

########## cavity calculation
cageMOL = "cage_ACIE_2006_45_901.mol2" #we need mol files as pdb file was not behaving OK with the aromatic

grid_spacing = 3
distance_threshold_for_90_deg_angle = 7
calculate_bfactor = True
compute_aromatic_contacts = False
compute_atom_contacts = True
distThreshold_atom_contacts = 5.0
dummy_atom_radii = 3
distanceFromCOMFactor = 0.0 #Distance from the center of mass factor 0-1 to consider spherical cavity, only useful for large cavities

threads_KDThree = 4

printLevel = 1 #print level 1 = normal print, print level 2 = print all info


######################################################################
###   Calculate cavity
######################################################################
print("\n--- Calculation of the cavity ---")

start_time = time.time()

## Load cage MOL file using RDKit
cagePDB = cageMOL.replace(".mol2", ".pdb")
cagePDBout1 = cagePDB.replace(".pdb", "_cavity.pdb")
rdkit_cage = rdkit.MolFromMol2File(cageMOL,removeHs=False)
#rdkit_cage = rdkit.MolFromPDBFile(cagePDB,sanitize=False,removeHs=False,proximityBonding=False)

#Calculate the center of mass of the cage using RDKit
atoms = [ atom for atom in rdkit_cage.GetAtoms() ]
numberOfAtoms = rdkit_cage.GetNumAtoms()
xyzPositionArray =  np.array([list( rdkit_cage.GetConformer().GetAtomPosition(i) ) for i in range(numberOfAtoms)])
center_of_mass = np.array(sum(atoms[i].GetMass()*xyzPositionArray[i] for i in range(numberOfAtoms)))/Descriptors.MolWt(rdkit_cage)
print("center_of_mass= ", center_of_mass)

# Assign the hydrophobic values to cage atoms from the library
listGlobalOut = list()
listGlobalOut2 = list() 
for i in hydrophValues:
    match = rdkit_cage.GetSubstructMatches(Chem.MolFromSmarts(hydrophValues[i][1]),False,False)
    listOut = list()
    for x in match:
        listOut.append(x[0]+1)
    listGlobalOut2.append(list(set(listOut)))
    if listOut:
        listGlobalOut.extend(list(set(listOut)))

numberOfAtoms = rdkit_cage.GetNumAtoms()

atomTypesList = []
for i in range(0,numberOfAtoms):
    atomTypesList.append([])

atomTypesListHydrophValues = []
for i in range(0,numberOfAtoms):
    atomTypesListHydrophValues.append([])

atomTypesMeanListHydrophValues = []

for i in range(0, len(listGlobalOut2)):
    if listGlobalOut2[i]:
        for j in listGlobalOut2[i]:          
            atomTypesList[j-1].append(i+1)

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
    if (printLevel == 2): print(atomSymbol, i+1, atomTypesList[i], valuesList, meanValuestList)



# Get cage xyx atom coordinates
xyzAtoms = np.array([rdkit_cage.GetConformer().GetAtomPosition(cage_atom.GetIdx()) for cage_atom in rdkit_cage.GetAtoms()], dtype=np.float64)
# Calculate neighbors with KDTree. Leaf_size affects the speed of a query and the memory required to store the constructed tree. 
# The amount of memory needed to store the tree scales as approximately n_samples / leaf_size.
kdtxyzAtoms = KDTree(xyzAtoms, leafsize = 20)

#We use the center of mass of the molecule as the pore center of mass
pore_center_of_mass = center_of_mass
distancesFromCOM = kdtxyzAtoms.query(center_of_mass, k = None, p = 2)
maxDistanceFromCOM = (distancesFromCOM[0][-1])
pore_radius = distancesFromCOM[0][0]*distanceFromCOMFactor

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

# Create a KDThree of the dummy atom gird
calculatedGirdKDTree = KDTree(calculatedGird.gridPosList, leafsize = 20)

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

#Create a KDTree for each atom type and store in a dict
KDTree_dict = {}
for atom_type in atom_type_list:
    KDTree_dict[atom_type] = KDTree(coords_dict[atom_type], leafsize = 20)

#Radius of the dummy D atoms in the cavity defined with a custom diameter
#vdwRdummy = rdkit.GetPeriodicTable().GetRvdw(1)
vdwRdummy = dummy_atom_radii
for i in calculatedGird.grid:
    dist1 = np.linalg.norm(np.array(pore_center_of_mass)-np.array(i.pos))
    if dist1 == 0:
        i.inside_cavity = 1
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
                    i.overlapping_with_cage = 1

                for atom_pos_index in xyzAtomsSet2[1]:
                    atom_pos = coords_dict[atom_type][atom_pos_index]
                    dist2 = np.linalg.norm(np.array(atom_pos)-np.array(i.pos))
                    vect2Norm = (np.array(atom_pos)-np.array(i.pos)) / dist2
                    angle = np.arccos(np.dot(vect1Norm, vect2Norm))
                    summAngles.append(angle)
                    distancesAngles.append(1/dist2)
        if (summAngles):
            averageSummAngles = np.average(summAngles, axis=None, weights=distancesAngles)
            averageSummAngles_deg = np.degrees(averageSummAngles)
            i.vector_angle = averageSummAngles_deg 
            if i.overlapping_with_cage == 0:
                if averageSummAngles_deg > 90:
                    i.inside_cavity = 1


# Create a KDThree of the dummy atoms
calculatedGirdContacts = []
for i in calculatedGird.grid:
    if i.inside_cavity == 1:
        calculatedGirdContacts.append(i.pos) 
calculatedGirdContactsKDTree = KDTree(calculatedGirdContacts, leafsize = 20)

# Calculate the numer of neighbors for each dummy atom. From 1 to 7 (as 1 contact is always, the atom itself)
for i, dummy_atom in enumerate(calculatedGird.grid):
    if dummy_atom.inside_cavity == 1:
        xyzDummySet = calculatedGirdContactsKDTree.query(dummy_atom.pos, k=None, p=2, distance_upper_bound=1.1*grid_spacing)
        dummy_atom.number_of_neighbors = len(xyzDummySet[1])

'''
##This part of the code can be removed as overlapping with the cage is calculated before
# Check the dummy atoms that overlap with the cage
for i in calculatedGird.grid:   
    for atom_type in atom_type_list:
        vdwR = vdwR_dict[atom_type][0]
        distThreshold = vdwR + vdwRdummy
        xyzAtomsSet2 = KDTree_dict[atom_type].query(i.pos, k = None, p = 2, distance_upper_bound = distThreshold, workers = threads_KDThree)
        if xyzAtomsSet2[0]:
            i.overlapping_with_cage = 1
        else:
            i.overlapping_with_cage = 0
'''

#Create an empty editable molecule in RDKit
rdkit_molRW = rdkit.RWMol()
rdkit_conf = rdkit.Conformer(1)

for i, dummy_atom in enumerate(calculatedGird.grid):
    if (dummy_atom.inside_cavity == 1 and dummy_atom.overlapping_with_cage == 0):
        # Remove isolated dummy atoms that only have one neighbors 
        if ( dummy_atom.number_of_neighbors > 1):
            if calculate_bfactor == True:
                bfactor_info = []                
                for atom_type in atom_type_list:
                    dist = KDTree_dict[atom_type].query(dummy_atom.pos, k = None, p = 2, distance_upper_bound = distThreshold_atom_contacts, workers = threads_KDThree)

                    if dist[0]:
                        for k in range (0, len(dist[0])):
                            if compute_atom_contacts == True:
                                bfactor_info.append(atomTypesMeanListHydrophValues[dist[1][k]]/(1+dist[0][k]))
                            elif compute_aromatic_contacts == True:
                                isAromatic = rdkit_cage.GetAtomWithIdx(atom_idx_dict[atom_type][dist[1][k]]).GetIsAromatic()
                                if isAromatic == True:
                                    bfactor_info.append(1/(1+dist[0][k]))
                                            
            #Add atom to the RDKit molecule
            indexNewAtom = rdkit_molRW.AddAtom(rdkit.Atom(1))            
            rdkit_conf.SetAtomPosition(indexNewAtom,Point3D(dummy_atom.x,dummy_atom.y,dummy_atom.z)) #Set xyz postion
            molecInfo = rdkit.AtomPDBResidueInfo()            
            molecInfo.SetName(' D') # the rdkit PDB name has incorrect whitespace
            molecInfo.SetResidueName('HOH')            
            molecInfo.SetResidueName(''.ljust(2-len('D'))+'HOH') # the rdkit PDB residue name has incorrect whitespace
            molecInfo.SetResidueNumber(indexNewAtom+1)
            molecInfo.SetIsHeteroAtom(False)
            molecInfo.SetOccupancy(1.0)
            if calculate_bfactor == True:
                molecInfo.SetTempFactor( sum(bfactor_info) )  #Set the b-factor value
            else:
                molecInfo.SetTempFactor( 0 )
            #molecInfo.SetTempFactor(dummy_atom.inside_cavity)
            #molecInfo.SetTempFactor(dummy_atom.d_from_center) # Save in the b-factor the distance from the center
            rdkit_molRW.GetAtomWithIdx(indexNewAtom).SetMonomerInfo(molecInfo)

rdkit_molRW.AddConformer(rdkit_conf, assignId=True)


print("Saving PDB file with the dummy cavity atoms")
rdkit.MolToPDBFile(rdkit_molRW, cagePDBout1)
#print("Saving MOL file with the dummy cavity atoms")
#rdkit.MolToMolFile(rdkit_molRW, cageMOLout1)

if len(rdkit_molRW.GetAtoms()) >0:
    ##Calculation of the volume of the cavity using RDKit
    #cageCavityVolume = rdkit.ComputeMolVolume(rdkit_molRW, confId=-1, gridSpacing=0.2, boxMargin=2.0)
    #print("Cage cavity volume (rdkit calc.) = ", cageCavityVolume, " A3")
    
    xyzPositionArrayMolRW =  np.array([list( rdkit_molRW.GetConformer().GetAtomPosition(i) ) for i in range(rdkit_molRW.GetNumAtoms())])
    
    cageCavityVolume = cavity_volume(xyzPositionArrayMolRW, radius=dummy_atom_radii, volume_grid_size=0.2)
    
    print("Cage cavity volume = ", cageCavityVolume, " A3")
    
    cavity_dummy_atoms = [ atom for atom in rdkit_molRW.GetAtoms() ]
    hydrophobicity_cavity_dummy_atoms = [ atom.GetMonomerInfo().GetTempFactor() for atom in cavity_dummy_atoms ]

    if (printLevel == 2): print("Individual hydrophobicity",hydrophobicity_cavity_dummy_atoms)
    average_cavity_hydrophobicity = np.mean(hydrophobicity_cavity_dummy_atoms)
    total_cavity_hydrophobicity = average_cavity_hydrophobicity*cageCavityVolume
    print("Average cavity hydrophobicity =", average_cavity_hydrophobicity,  " A^-3")
    print("Total cavity hydrophobicity =", total_cavity_hydrophobicity )   
  
    cage_name = cagePDB.replace(".pdb", "")
    file1 = open(cage_name+"_cavity_vol_calc.txt","w+")
    file1.write("Cage cavity volume = " + str(cageCavityVolume) + " A3\n")
    file1.write("Average cavity hydrophobicity = " + str(total_cavity_hydrophobicity) + "\n")
    file1.write("Total cavity hydrophobicity = " + str(total_cavity_hydrophobicity) + " A^-3\n")
    file1.close()
else:
    print("Cavity with no volume")



print("--- Total time %s seconds ---" % (time.time() - start_time))

######################################################################
###        Open the saved cavity in PyMol
######################################################################

# pymol launching: quiet (-q)
import pymol
pymol.pymol_argv = ['pymol','-q']
pymol.finish_launching()
cmd = pymol.cmd
cmd.load(cagePDB, "cage")
cmd.set('valence', 0)
cmd.load(cagePDBout1, "cavity")

cmd.alter('name D', 'vdw="' + str(dummy_atom_radii) + '"')

#cmd.show_as("nonbonded", selection="cavity")
#cmd.show_as("surface", selection="cavity")
cmd.show_as("spheres", selection="cavity")
cmd.show_as("spheres", selection="cage")

cmd.spectrum("b", selection="cavity")
#cmd.set("surface_color", "withe", selection="cavity")
#cmd.set("transparency", 0.5, cage_name)

cmd.ramp_new("ramp", "cavity", [min(hydrophobicity_cavity_dummy_atoms),max(hydrophobicity_cavity_dummy_atoms)], "rainbow")


cmd.clip("atoms", 5, "All")
cmd.orient("cage")
cmd.zoom("cage")

cmd.save(cagePDBout1.replace(".pdb", ".pse"))
