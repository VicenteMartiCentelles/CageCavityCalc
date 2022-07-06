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

from hydrophobicity_values import hydrophValuesGhose1998, hydrophValuesCrippen1999
from calculations import sum_grid_volume, assignHydrophobicValuesToCageAtoms
from cavity_classes import GridPoint, CageGrid


from input_output import read_positions_and_atom_names_from_file, print_to_file, read_cgbind, read_mdanalysis
#from old.cavity_calculator_v08 import cagePDB


class cavity():
    def __init__(self):
        '''
        Initilalize the methodology, set up all the variables
        '''

        self.grid_spacing = 1 # spacing between grid points
        self.distance_threshold_for_90_deg_angle = 7
        calculate_bfactor = True
        compute_aromatic_contacts = False
        compute_hydrophobicity = True
        distance_function = "Fauchere" #Audry, Fauchere, Fauchere2 // Distance function of the hydrophobic potential, OnlyValues to do not consider the distance (only used for testing purpouses)
        distThreshold_atom_contacts = 20.0
        self.dummy_atom_radii = 1 #radius of the dummy atom use for the cavity calculations
        self.distanceFromCOMFactor = 0.4 #Distance from the center of mass factor 0-1 to consider spherical cavity, only useful for large cavities

        self.threads_KDThree = 4

        printLevel = 1 #print level 1 = normal print, print level 2 = print all info

        #hydrophValues = hydrophValuesGhose1998
        hydrophValues = hydrophValuesCrippen1999


        # proporties of the loaded cage:
        self.positions = None
        self.atom_names = None
        self.atom_masses = None
        self.atom_vdw = None
        self.n_atoms = 0

        self.dummy_atoms_positions = []
        self.dummy_atoms_temperature = None

        #self.hydrophobicity = hydrophobicity

        self.volume = None

    def read_file(self, filename):
        self.positions, self.atom_names, self.atom_masses, self.atom_vdw = read_positions_and_atom_names_from_file(filename)
        self.n_atoms = len(self.positions)
        # Reset volume if the file is read  (in case it was calculated for the diffrent file)
        self.volume = None
        self.dummy_atoms = []

        ########## cavity calculation
        '''
        cageMOL = "test.mol2"  # we need mol files as pdb file was not behaving OK with the aromatic
        os.chdir("./test")
        print("Current Working Directory ", os.getcwd())
        
        ## Load cage MOL file using RDKit
        cagePDB = cageMOL.replace(".mol2", ".pdb")
        cagePDBout1 = cagePDB.replace(".pdb", "_cavity.pdb")
        rdkit_cage = rdkit.MolFromMol2File(cageMOL, removeHs=False)
        '''

    def read_cgbind(self, cgbind_cage):
        self.positions, self.atom_names, self.atom_masses, self.atom_vdw = read_cgbind(cgbind_cage)
        self.n_atoms = len(self.positions)

        # Reset volume if the file is read  (in case it was calculated for the diffrent file)
        self.volume = None
        self.dummy_atoms = []


    def read_mdanalysis(self, syst):
        self.positions, self.atom_names, self.atom_masses, self.atom_vdw = read_mdanalysis(syst)
        self.n_atoms = len(self.positions)

        # Reset volume if the file is read  (in case it was calculated for the diffrent file)
        self.volume = None
        self.dummy_atoms = []


    def volume(self):
        if self.volume is None:
            self.calculate_volume()
        else:
            return self.volume

    def calculate_volume(self):

        print("\n--- Calculation of the cavity ---")

        # TODO check the volume and check if positions ant atoms names not empty

        start_time = time.time()

        pore_center_of_mass, pore_radius = self.calculate_center_and_radius()

        calculatedGird = self.set_up_grid(pore_center_of_mass)
        #self.create_empty_molecule()

        self.find_dummies_inside_cavity(calculatedGird, pore_center_of_mass, pore_radius)
        self.volume = self.sum_up_volume()
        print(f"--- Total time {(time.time() - start_time):.0f} seconds ---" )

    def hydrophobicity(self):
        None #TODO

    def pymol(self):
        None #TODO

    def calculate_center_and_radius(self):
        #Calculate the center of mass of the cage using RDKit
        pore_center_of_mass = np.array(sum(self.atom_masses[i]*self.positions[i] for i in range(self.n_atoms)))/sum(self.atom_masses)
        print("center_of_mass= ", pore_center_of_mass)

        #We use the center of mass of the molecule as the pore center of mass
        kdtxyzAtoms = KDTree(self.positions, leafsize = 20) #Calculate the KDTree of the cage atom positions
        distancesFromCOM = kdtxyzAtoms.query(pore_center_of_mass, k = None, p = 2)
        maxDistanceFromCOM = (distancesFromCOM[0][-1])
        pore_radius = distancesFromCOM[0][0]*self.distanceFromCOMFactor

        return pore_center_of_mass, pore_radius

        # Assign the hydrophobic values to cage atoms from the library
        #HYDRO atomTypesMeanListHydrophValues, atomTypesValuesListHydrophValues, atomTypesInfoAtomSymbol, atomTypesInfoAtomGlobalIndex, atomTypesAssignemet = assignHydrophobicValuesToCageAtoms(rdkit_cage, hydrophValues, printLevel)

        '''HYDRO
                atom_symbols_in_molecule = []
                for atom_symbol in atomTypesInfoAtomSymbol:
                    if atom_symbol not in atom_symbols_in_molecule:
                        atom_symbols_in_molecule.append(atom_symbol)
                
                atomTypes_HydrophValues_dict = {}
                for j, atom_type in enumerate(atom_symbols_in_molecule):
                    k = 1
                    atomTypes_dict = {}
                    for i, atom_symbol in enumerate(atomTypesInfoAtomSymbol):
                        if atom_type == atom_symbol: 
                            atomTypes_dict[k] = i + 1
                            k = k + 1
                    atomTypes_HydrophValues_dict[atom_type] = atomTypes_dict
        
                print(atomTypes_HydrophValues_dict)
           
        if (printLevel == 2): 
            fileExtraData = open(cagePDB.replace(".pdb", "_extra_data.txt"),"w+")
            fileExtraData.write("Asignement of hydrophobic values to cage atoms from tables\n")
            for i, atom_symbol in enumerate(atomTypesInfoAtomSymbol):
                fileExtraData.write( str(atom_symbol) + str(i+1)+ " Atom type:" + str(atomTypesAssignemet[i]) + str(atomTypesValuesListHydrophValues[i]) + str(atomTypesMeanListHydrophValues[i]) + "\n")
            fileExtraData.write("*****************\n\n")
        '''

    def set_up_grid(self, pore_center_of_mass):
        # Calculate the cage bounding box using the x,y,z axis
        #box_size = maxDistanceFromCOM ## Not used, in this previous version we used a square box
        x_coordinates, y_coordinates, z_coordinates = zip(*self.positions)
        box_size = [min(x_coordinates),max(x_coordinates),min(y_coordinates),max(y_coordinates),min(z_coordinates),max(z_coordinates)]
        print("Box x min/max= ", box_size[0], box_size[1])
        print("Box y min/max= ", box_size[2], box_size[3])
        print("Box z min/max= ", box_size[4], box_size[5])

        calculatedGird = CageGrid(pore_center_of_mass, box_size, delta = 0, grid_spacing = self.grid_spacing)
        return calculatedGird

    def create_dummy_molecule(self):
        if(printLevel == 2):
            #Create an empty editable molecule in RDKit, to save the dummy atom box
            rdkit_molRW = rdkit.RWMol()
            rdkit_conf = rdkit.Conformer(1)

            for i, dummy_atom in enumerate(calculatedGird.grid):
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
                molecInfo.SetTempFactor(1.0)
                rdkit_molRW.GetAtomWithIdx(indexNewAtom).SetMonomerInfo(molecInfo)

            rdkit_molRW.AddConformer(rdkit_conf, assignId=True)

            print("Saving PDB file with the box of dummy cavity atoms")
            rdkit.MolToPDBFile(rdkit_molRW, cagePDB.replace(".pdb", "_box_cavity.pdb"))

    def find_dummies_inside_cavity(self, calculatedGird, pore_center_of_mass, pore_radius):
        #Create a KDTree for each atom type and store in a dict
        atom_type_list = []
        coords_dict = {}
        vdwR_dict = {}
        atom_idx_dict = {}

        for atom_idx, cage_name in enumerate(self.atom_names):
            atom_type = cage_name
            pos = self.positions[atom_idx]

            if atom_type not in atom_type_list:
                atom_type_list.append(atom_type)

                vdwR_dict.setdefault(atom_type, []).append(self.atom_vdw[atom_idx])

            coords_dict.setdefault(atom_type, []).append(list(pos))
            atom_idx_dict.setdefault(atom_type, []).append(atom_idx)

        KDTree_dict = {}
        for atom_type in atom_type_list:
            KDTree_dict[atom_type] = KDTree(coords_dict[atom_type], leafsize = 20)

        #Calculate the dummy atoms that form the cavity using the angle method, also check that atoms do not overlap with cage
        vdwRdummy = self.dummy_atom_radii #Radius of the dummy D atoms in the cavity defined with a custom diameter
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
                    xyzAtomsSet2 = KDTree_dict[atom_type].query(i.pos, k = None, p = 2,
                                                                distance_upper_bound = self.distance_threshold_for_90_deg_angle)
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
                            distancesAngles.append(1/(1+dist2))
                if (summAngles):
                    averageSummAngles = np.average(summAngles, axis=None, weights=distancesAngles)
                    averageSummAngles_deg = np.degrees(averageSummAngles)
                    i.vector_angle = averageSummAngles_deg
                    if i.overlapping_with_cage == 0:
                        if averageSummAngles_deg > 90:
                            i.inside_cavity = 1

        # Create a KDThree of the cavity dummy atoms
        calculatedGirdContacts = []
        for i in calculatedGird.grid:
            if i.inside_cavity == 1:
                calculatedGirdContacts.append(i.pos)
        calculatedGirdContactsKDTree = KDTree(calculatedGirdContacts, leafsize = 20)

        # Calculate the number of neighbors for each dummy atom. From 1 to 7 (as 1 contact is always, the atom itself)
        for i, dummy_atom in enumerate(calculatedGird.grid):
            if dummy_atom.inside_cavity == 1:
                xyzDummySet = calculatedGirdContactsKDTree.query(dummy_atom.pos, k=None, p=2, distance_upper_bound=1.1*self.grid_spacing)
                dummy_atom.number_of_neighbors = len(xyzDummySet[1])


        #Create an empty editable molecule in RDKit, store bfactor data then save as PDB
        '''
        rdkit_molRW = rdkit.RWMol()
        rdkit_conf = rdkit.Conformer(1)
        '''

        for i, dummy_atom in enumerate(calculatedGird.grid):
            if (dummy_atom.inside_cavity == 1 and dummy_atom.overlapping_with_cage == 0):
                # Remove isolated dummy atoms that only have one neighbors
                if ( dummy_atom.number_of_neighbors > 1):

                    # This is some extra stuff
                    '''
                    if calculate_bfactor == True:
                        bfactor_info = []
                        if (printLevel == 2): fileExtraData.write( "-------------\n" )
                        for atom_type in atom_type_list:
                            dist = KDTree_dict[atom_type].query(dummy_atom.pos, k = None, p = 2, distance_upper_bound = distThreshold_atom_contacts, workers = threads_KDThree)
                            if (printLevel == 2):
                                fileExtraData.write( str(atom_type) + ": " + str(dist[0]) + "\n" )
                                fileExtraData.write( str(atom_type) + " Index: "+ str([x+1 for x in dist[1]]) + "\n" )# Add 1 to key as atom number starts with 1 and list number with 0
                                fileExtraData.write( "Global Index: " + str([atomTypes_HydrophValues_dict[atom_type].get(key+1) for key in dist[1]]) + "\n" ) # Add 1 to key as atom number starts with 1 and list number with 0

                            if dist[0]:
                                for k in range (0, len(dist[0])):
                                    Hydroph_Value = atomTypesMeanListHydrophValues[atomTypes_HydrophValues_dict[atom_type][1+dist[1][k]] -1 ] # we need to add +1 to the k index as atom numbering is starts at 1 and lists index at 0. After taht, we need to add -1 to the atom index from the dict to obtain the value from the list that starts from 0
                                    if compute_hydrophobicity == True:
                                        if distance_function == "Audry":
                                            bfactor_info.append(Hydroph_Value/(1+dist[0][k]))
                                        elif distance_function == "Fauchere":
                                            bfactor_info.append(Hydroph_Value*np.exp(-1*dist[0][k]))
                                        elif distance_function == "Fauchere2":
                                            bfactor_info.append(Hydroph_Value*np.exp(-1/2*dist[0][k]))
                                        elif distance_function == "OnlyValues":
                                            bfactor_info.append( Hydroph_Value )
                                    elif compute_aromatic_contacts == True:
                                        isAromatic = rdkit_cage.GetAtomWithIdx(atom_idx_dict[atom_type][dist[1][k]]).GetIsAromatic()
                                        if isAromatic == True:
                                            bfactor_info.append(1/(1+dist[0][k]))
                    '''

                    #Add atom to the RDKit molecule

                    self.dummy_atoms_positions.append([dummy_atom.x,dummy_atom.y,dummy_atom.z])
                    '''
                    indexNewAtom = rdkit_molRW.AddAtom(rdkit.Atom(1))
                    rdkit_conf.SetAtomPosition(indexNewAtom,Point3D(dummy_atom.x,dummy_atom.y,dummy_atom.z)) #Set xyz postion
                    molecInfo = rdkit.AtomPDBResidueInfo()
                    molecInfo.SetName(' D') # the rdkit PDB name has incorrect whitespace
                    molecInfo.SetResidueName('HOH')
                    molecInfo.SetResidueName(''.ljust(2-len('D'))+'HOH') # the rdkit PDB residue name has incorrect whitespace
                    molecInfo.SetResidueNumber(indexNewAtom+1)
                    molecInfo.SetIsHeteroAtom(False)
                    molecInfo.SetOccupancy(1.0)
                    '''
                    '''
                    if calculate_bfactor == True:
                        molecInfo.SetTempFactor( sum(bfactor_info) )  #Set the b-factor value
                        if (printLevel == 2):
                            fileExtraData.write( "Dummy atom index = " + str(indexNewAtom+1) + "\n" )
                            fileExtraData.write( "Sum Hydrop. =" + str(sum(bfactor_info)) + " Individual = " + str(bfactor_info)+ "\n" )
                    else:
                        molecInfo.SetTempFactor( 0 )
                    '''

                    #molecInfo.SetTempFactor(dummy_atom.inside_cavity)
                    #molecInfo.SetTempFactor(dummy_atom.d_from_center) # Save in the b-factor the distance from the center
                    '''
                    rdkit_molRW.GetAtomWithIdx(indexNewAtom).SetMonomerInfo(molecInfo)
                    '''

        '''
        rdkit_molRW.AddConformer(rdkit_conf, assignId=True)

        if (printLevel == 2): fileExtraData.close()

        print("Saving PDB file with the dummy cavity atoms")
        rdkit.MolToPDBFile(rdkit_molRW, cagePDBout1)
        #print("Saving MOL file with the dummy cavity atoms")
        #rdkit.MolToMolFile(rdkit_molRW, cageMOLout1)
        '''

    def sum_up_volume(self):
        #Calculate the volume of the cavity, print and save
        if len(self.dummy_atoms_positions) >0:
            cageCavityVolume = sum_grid_volume(np.array(self.dummy_atoms_positions),
                                               radius=self.dummy_atom_radii,
                                               volume_grid_size=self.grid_spacing/2)

            print("Cage cavity volume = ", cageCavityVolume, " A3")
            return cageCavityVolume

            '''
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
            file1.write("Average cavity hydrophobicity = " + str(average_cavity_hydrophobicity) + " A^-3\n")
            file1.write("Total cavity hydrophobicity = " + str(total_cavity_hydrophobicity) + "\n")
            file1.close()
            '''

        else:
            print("Cavity with no volume")

    #def to_pdb_file(self):

    def print_to_file(self, filename): #TODO check if cavity
        if len(self.dummy_atoms_positions)==0:
            print("No cavity, saving just the input!") #TODO
            print_to_file(filename, self.positions, self.atom_names)
        else:
            positions = np.vstack([self.positions, self.dummy_atoms_positions])
            #positions2 = np.hstack([self.positions, self.dummy_atoms_positions])

            atom_names = np.append(self.atom_names, np.array(['D']*len(self.dummy_atoms_positions)))
            print_to_file(filename, positions, atom_names)





