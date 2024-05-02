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
from scipy.spatial import distance_matrix

from sklearn.cluster import DBSCAN


from CageCavityCalc.data import hydrophValuesGhose1998, hydrophValuesCrippen1999, vdw_radii
from CageCavityCalc.calculations import sum_grid_volume
from CageCavityCalc.grid_classes import GridPoint, CageGrid
from CageCavityCalc.hydrophobicity import assignHydrophobicValuesToCageAtoms, calc_single_hydrophobicity
from CageCavityCalc.electrostatics import calculate_partial_charges
from CageCavityCalc.input_output import read_positions_and_atom_names_from_file, read_positions_and_atom_names_from_array, print_to_file, read_cgbind, read_mdanalysis, print_pymol_file, convert_to_mol2
from CageCavityCalc.window_size import get_max_escape_sphere
from CageCavityCalc.log import logger

class cavity():
    def __init__(self):
        '''
        Initilalize the methodology, set up all the variables
        '''

        self.grid_spacing = 1 # spacing between grid points

        self.distance_threshold_for_90_deg_angle = 5
        calculate_bfactor = True
        compute_aromatic_contacts = False
        compute_hydrophobicity = True
        self.dummy_atom_radii = 1 #radius of the dummy atom use for the cavity calculations
        self.clustering_to_remove_cavity_noise = "false" #"Size" select largest cluster, "Dist" select closte cluster to the center
        
        # With the modified code we can always use distanceFromCOMFactor = 1
        self.distanceFromCOMFactor = 1 #Distance from the center of mass factor 0-1 to consider spherical cavity, only useful for large cavities

        self.threads_KDThree = 4

        printLevel = 1 #print level 1 = normal print, print level 2 = print all info
        
        self.hydrophMethod = "Ghose" # Hydrophobic calculation method "Ghose" or "Crippen"


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

        self.filename = None # original file of the cage (might be needed to conversion to rdkit in case of hydrophobicity calculation)
        self.compute_hydrophobicity = False
        self.KDTree_dict = None
        self.distThreshold_atom_contacts = 20.0
        self.distance_function = "Fauchere" #Audry, Fauchere, Fauchere2 // Distance function of the hydrophobic potential, OnlyValues to do not consider the distance (only used for testing purpouses)
        # Calculated properties (since they are easy to calculate
        self.hydrophobicity = []
        self.aromatic_constacts = []
        self.solvent_accessibility = []
        self.esp_grid = []


    def read_file(self, filename):
        logger.info("--- Reading the file ---")

        self.filename = filename
        self.positions, self.atom_names, self.atom_masses, self.atom_vdw = read_positions_and_atom_names_from_file(filename)
        self.n_atoms = len(self.positions)
        # Reset volume if the file is read  (in case it was calculated for the diffrent file)
        self.volume = None
        self.dummy_atoms = []

    def read_pos_name_array(self, positions, names):
        logger.info("--- Reading the arrays ---")
        self.positions, self.atom_names, self.atom_masses, self.atom_vdw = read_positions_and_atom_names_from_array(np.array(positions), np.array(names))
        self.n_atoms = len(self.positions)
        # Reset volume if the file is read  (in case it was calculated for the diffrent file)
        self.volume = None
        self.dummy_atoms = []

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

        logger.info("--- Calculation of the cavity ---")

        assert self.n_atoms > 0, "no atoms"

        start_time = time.time()

        pore_center_of_mass, pore_radius = self.calculate_center_and_radius()

        calculatedGird = self.set_up_grid(pore_center_of_mass)
        #self.create_empty_molecule()

        self.find_dummies_inside_cavity(calculatedGird, pore_center_of_mass, pore_radius)
        self.volume = self.sum_up_volume()
        logger.info(f"--- Total time {(time.time() - start_time):.0f} seconds ---" )
        return self.volume

    def calculate_center_of_mass(self):
        pore_center_of_mass = np.array(sum(self.atom_masses[i]*self.positions[i] for i in range(self.n_atoms)))/sum(self.atom_masses)
        return pore_center_of_mass.tolist()
        
    def calculate_center_and_radius(self):
        #Calculate the center of mass of the cage using RDKit
        pore_center_of_mass = np.array(sum(self.atom_masses[i]*self.positions[i] for i in range(self.n_atoms)))/sum(self.atom_masses)
        logger.info(f"center_of_mass= {pore_center_of_mass:}")

        #We use the center of mass of the molecule as the pore center of mass
        kdtxyzAtoms = KDTree(self.positions, leafsize = 20) #Calculate the KDTree of the cage atom positions
        distancesFromCOM = kdtxyzAtoms.query(pore_center_of_mass, k = self.n_atoms, p = 2)
        #print("++++++++++++++++++++++++++++++++++++++++++")
        #print("distancesFromCOM: ")
        #print(distancesFromCOM)
        #print(len(distancesFromCOM))
        maxDistanceFromCOM = (distancesFromCOM[0][-1])
        pore_radius = distancesFromCOM[0][0]*self.distanceFromCOMFactor
        '''
        #Determine the atom types and correct the pore radius  with the atom radius of the cage. We use the max radius of the different cage atom types
        atomTypesInCage = list(dict.fromkeys(self.atom_names))
        print(atomTypesInCage)
        atomTypesInCageRadius=[*map(vdw_radii.get, atomTypesInCage)]
        print(atomTypesInCageRadius)
        atomTypesInCageRadius = list(filter(lambda item: item is not None, atomTypesInCageRadius))
        print(max(atomTypesInCageRadius))
        pore_radius = pore_radius - max(atomTypesInCageRadius) - self.dummy_atom_radii
        '''
        #Determine the atom type of the calculated distance and correct the pore radius  with the atom radius and the dummy atom radius
        pore_radius = pore_radius - 1.01*vdw_radii[self.atom_names[distancesFromCOM[1][0]]] - 1.01*self.dummy_atom_radii
        if pore_radius < 0:
            pore_radius = 0

        return pore_center_of_mass, pore_radius

        # Assign the hydrophobic values to cage atoms from the library

        '''HYDRO

           
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
        logger.info(f"Box x min/max= {box_size[0]:f}, {box_size[1]:f}")
        logger.info(f"Box y min/max= {box_size[2]:f}, {box_size[3]:f}")
        logger.info(f"Box z min/max= {box_size[4]:f}, {box_size[5]:f}")

        calculatedGird = CageGrid(pore_center_of_mass, box_size, delta = 0, grid_spacing = self.grid_spacing)
        return calculatedGird

    def find_dummies_inside_cavity(self, calculatedGird, pore_center_of_mass, pore_radius):
        #Create a KDTree for each atom type and store in a dict
        self.atom_type_list = []
        coords_dict = {}
        vdwR_dict = {}
        self.atom_idx_dict = {}

        for atom_idx, cage_name in enumerate(self.atom_names):
            atom_type = cage_name
            pos = self.positions[atom_idx]

            if atom_type not in self.atom_type_list:
                self.atom_type_list.append(atom_type)

                vdwR_dict.setdefault(atom_type, []).append(self.atom_vdw[atom_idx])

            coords_dict.setdefault(atom_type, []).append(list(pos))
            self.atom_idx_dict.setdefault(atom_type, []).append(atom_idx)

        self.KDTree_dict = {}
        for atom_type in self.atom_type_list:
            self.KDTree_dict[atom_type] = KDTree(coords_dict[atom_type], leafsize = 20)



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
                for atom_type in self.atom_type_list:
                    xyzAtomsSet2 = self.KDTree_dict[atom_type].query(i.pos, k = self.n_atoms, p = 2,
                                                                distance_upper_bound = self.distance_threshold_for_90_deg_angle)




                    #print("xyzAtomsSet2[0]:",xyzAtomsSet2[0])
                    for index_atom_iter in range(0, len(xyzAtomsSet2[0])):
                        if xyzAtomsSet2[0][index_atom_iter] == math.inf:
                            index_inf = index_atom_iter
                            xyzAtomsSet3_dist = xyzAtomsSet2[0][:index_inf]
                            xyzAtomsSet3_index = xyzAtomsSet2[1][:index_inf]
                            break

                    #print("xyzAtomsSet3_dist:",xyzAtomsSet3_dist)
                    if len(xyzAtomsSet3_dist) > 0:
                        #print("if xyzAtomsSet2[1]:",xyzAtomsSet2[0])
                        vdwR = vdwR_dict[atom_type][0]
                        distThreshold = vdwR + vdwRdummy
                        if xyzAtomsSet3_dist[0] < distThreshold:
                            i.overlapping_with_cage = 1

                        for atom_pos_index in xyzAtomsSet3_index:
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
                xyzDummySet = calculatedGirdContactsKDTree.query_ball_point(dummy_atom.pos, r=1.1*self.grid_spacing, p=2)
                dummy_atom.number_of_neighbors = len(xyzDummySet)
        
        if self.clustering_to_remove_cavity_noise != "false":
            cavity_dummy_atoms_positions = []
            cavity_dummy_atoms_index = []
            for i, dummy_atom in enumerate(calculatedGird.grid):
                if dummy_atom.inside_cavity == 1:
                    cavity_dummy_atoms_positions.append([dummy_atom.x,dummy_atom.y,dummy_atom.z])
                    cavity_dummy_atoms_index.append(i)
        
            clusters = DBSCAN(eps=self.grid_spacing*1.1).fit(np.array(cavity_dummy_atoms_positions))
            cluster_labels = clusters.labels_
            number_of_clusters = len( np.unique(cluster_labels) )
            number_of_noise = np.sum(np.array(cluster_labels) == -1, axis=0)
            largest_cluster = max(set(cluster_labels.tolist()), key = cluster_labels.tolist().count)
            print('Clusters labels:', np.unique(cluster_labels))
            print('Number of clusters: %d' % number_of_clusters)
            print('Number of noise points: %d' % number_of_noise)

            #Create a dict with the calculatedGird index and cluster label
            cavity_dummy_atoms_clusters = {cavity_dummy_atoms_index[i]: cluster_labels[i] for i in range(len(cavity_dummy_atoms_index))}
        

        if self.clustering_to_remove_cavity_noise == "dist":
            inv_map_cavity_dummy_atoms_clusters = {}
            for k, v in cavity_dummy_atoms_clusters.items():
                inv_map_cavity_dummy_atoms_clusters[v] = inv_map_cavity_dummy_atoms_clusters.get(v, []) + [k]
            
            all_grid_dummy_atoms_positions = []
            for i, dummy_atom in enumerate(calculatedGird.grid):
                all_grid_dummy_atoms_positions.append([dummy_atom.x,dummy_atom.y,dummy_atom.z])
            all_grid_dummy_atoms_positions = np.array(all_grid_dummy_atoms_positions)
            
            cluster_centroids = []
            cluster_centroids_distance_to_COM = []
            center_of_mass = self.calculate_center_of_mass()
            for i in np.unique(cluster_labels):
                cluster_centroids.append(np.mean(all_grid_dummy_atoms_positions[inv_map_cavity_dummy_atoms_clusters[i]], axis=0))
            
            for i in cluster_centroids:
                cluster_centroids_distance_to_COM.append(np.linalg.norm(i-center_of_mass))
            
            print(cluster_centroids_distance_to_COM)
            min_cluster_centroids_distance_to_COM = cluster_centroids_distance_to_COM[1:].index(min(cluster_centroids_distance_to_COM[1:]))

        if self.clustering_to_remove_cavity_noise == "size":
            print('Saving only largest cluster of the cavity')
            print('Number of points in the cluster: %d' % cluster_labels.tolist().count(largest_cluster) )
        elif self.clustering_to_remove_cavity_noise == "dist":
            print('Saving only cluster closest to the center of mass of the cavity')
            print('Number of points in the cluster: %d' % cluster_labels.tolist().count(min_cluster_centroids_distance_to_COM) )
            
        for i, dummy_atom in enumerate(calculatedGird.grid):
            if (dummy_atom.inside_cavity == 1 and dummy_atom.overlapping_with_cage == 0):
                # Remove isolated dummy atoms that only have one neighbors
                if ( (self.clustering_to_remove_cavity_noise == "false" and dummy_atom.number_of_neighbors > 1) or (self.clustering_to_remove_cavity_noise == "size" and dummy_atom.number_of_neighbors > 1 and cavity_dummy_atoms_clusters[i] == largest_cluster) or (self.clustering_to_remove_cavity_noise == "dist" and dummy_atom.number_of_neighbors > 1 and cavity_dummy_atoms_clusters[i] == min_cluster_centroids_distance_to_COM)):
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

                    '''

                    #Add atom to the RDKit molecule

                    self.dummy_atoms_positions.append([dummy_atom.x,dummy_atom.y,dummy_atom.z])

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

        if (printLevel == 2): fileExtraData.close()

        print("Saving PDB file with the dummy cavity atoms")
        rdkit.MolToPDBFile(rdkit_molRW, cagePDBout1)
        #print("Saving MOL file with the dummy cavity atoms")
        #rdkit.MolToMolFile(rdkit_molRW, cageMOLout1)
        '''

    def sum_up_volume(self):
        logger.info("Summing the volume")

        #Calculate the volume of the cavity, print and save
        if len(self.dummy_atoms_positions) >0:
            cageCavityVolume = sum_grid_volume(np.array(self.dummy_atoms_positions),
                                               radius=self.dummy_atom_radii,
                                               volume_grid_size=self.grid_spacing/2)

            logger.info(f"Cage cavity volume = {cageCavityVolume:.2f}  A3")
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
            logger.info("Cavity with no volume")


    def get_property_values(self, property_name):
        property_values = None

        if property_name == "aromaticity" or property_name == "aromatic_contast" or property_name == "a":
            property_values = np.append(np.zeros((len(self.positions))), self.aromatic_constacts)
        elif property_name == "solvent_accessibility" or property_name == "solvent" or property_name == "s":
            property_values = np.append(np.zeros((len(self.positions))), self.solvent_accessibility)
        elif property_name == "hydrophobicity" or property_name == "hydro" or property_name == "h":
            property_values = np.append(np.zeros((len(self.positions))), self.hydrophobicity)
        elif property_name == "electrostatics" or property_name == "esp":
            property_values = np.append(np.zeros((len(self.positions))), self.esp_grid)

        return property_values

    def print_to_file(self, filename, property_name = None):
        logger.info("Printing to file")

        assert self.n_atoms > 0, "no atoms"

        if len(self.dummy_atoms_positions)==0:
            logger.info("No cavity, saving just the input!")
            print_to_file(filename, self.positions, self.atom_names)
        else:
            property_values = self.get_property_values(property_name)

            positions = np.vstack([self.positions, self.dummy_atoms_positions])

            atom_names = np.append(self.atom_names, np.array(['D']*len(self.dummy_atoms_positions)))
            print_to_file(filename, positions, atom_names, property_values)


    def print_to_pymol(self, filename, property_name = None):

        logger.info("Printing to pymol")


        property_values = self.get_property_values(property_name)

        #Firstly we need to save to pdb file
        positions = np.vstack([self.positions, self.dummy_atoms_positions])
        atom_names = np.append(self.atom_names, np.array(['D']*len(self.dummy_atoms_positions)))


        print_to_file(filename[:filename.find('.')]+".pdb", positions, atom_names, property_values)
        #Than pymol:
        if property_values is not None:
            print_pymol_file(filename, property_values[len(self.positions):], self.dummy_atom_radii)
        else:
            print_pymol_file(filename, None, self.dummy_atom_radii)

    def calculate_hydrophobicity(self):

        logger.info("--- Calculation of the hydrophobicity ---")

        if self.hydrophMethod == "Ghose":
            self.hydrophValues = hydrophValuesGhose1998
        elif self.hydrophMethod == "Crippen":
            self.hydrophValues = hydrophValuesCrippen1999

        assert self.volume is not None, "Cavity not calculated (use calculate_volume())" # Not sure, maybe just check and calculate volume (?)

        self.compute_hydrophobicity = True

        if self.filename is not None and self.filename.endswith('.mol2'):
            rdkit_cage = rdkit.MolFromMol2File(self.filename, removeHs=False)
        else:
            # use openbabel:
            logger.info("No .mol2 file, we will convert it using openbabel")
            rdkit_cage = convert_to_mol2(self.positions, self.atom_names)



        atomTypesMeanListHydrophValues, atomTypesValuesListHydrophValues, atomTypesInfoAtomSymbol, atomTypesInfoAtomGlobalIndex, atomTypesAssignemet = assignHydrophobicValuesToCageAtoms(rdkit_cage, self.hydrophValues)
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

        for dummy_atom in self.dummy_atoms_positions:
            #if (printLevel == 2): fileExtraData.write("-------------\n")
            dummy_hydrophobicity = []
            dummy_aromatic_contacts = []
            dummy_solvent_accessibility = []

            for atom_type in self.atom_type_list:
                dist = self.KDTree_dict[atom_type].query(dummy_atom, k=self.n_atoms, p=2,
                                                    distance_upper_bound=self.distThreshold_atom_contacts)
                '''
                if (printLevel == 2):
                    fileExtraData.write(str(atom_type) + ": " + str(dist[0]) + "\n")
                    fileExtraData.write(str(atom_type) + " Index: " + str(
                        [x + 1 for x in dist[1]]) + "\n")  # Add 1 to key as atom number starts with 1 and list number with 0
                    fileExtraData.write("Global Index: " + str([atomTypes_HydrophValues_dict[atom_type].get(key + 1) for key in
                                                                dist[1]]) + "\n")  # Add 1 to key as atom number starts with 1 and list number with 0
                '''
                for index_atom_iter in range(0, len(dist[0])):
                    if dist[0][index_atom_iter] == math.inf:
                        index_inf = index_atom_iter
                        dist_dist = dist[0][:index_inf]
                        dist_index = dist[1][:index_inf]
                        break

                if len(dist_dist) > 0:
                    for k in range(0, len(dist_dist)):

                        # calculate hydrophobicity
                        Hydroph_Value = atomTypesMeanListHydrophValues[atomTypes_HydrophValues_dict[atom_type][1 + dist_index[
                            k]] - 1]  # we need to add +1 to the k index as atom numbering is starts at 1 and lists index at 0. After taht, we need to add -1 to the atom index from the dict to obtain the value from the list that starts from 0

                        hydro = calc_single_hydrophobicity(dist_dist[k], Hydroph_Value, self.distance_function)
                        dummy_hydrophobicity.append(hydro)

                        #calculate aromatic contact
                        isAromatic = rdkit_cage.GetAtomWithIdx(self.atom_idx_dict[atom_type][dist_index[k]]).GetIsAromatic()
                        if isAromatic == True:
                            dummy_aromatic_contacts.append(1 / (1 + dist_dist[k]))
                        else:
                            dummy_aromatic_contacts.append(0)

                        #calculate solvent accessibility
                        dummy_solvent_accessibility.append(1/dist_dist[k])

            self.hydrophobicity.append(np.sum(dummy_hydrophobicity))
            self.aromatic_constacts.append(np.sum(dummy_aromatic_contacts))
            self.solvent_accessibility.append(np.sum(dummy_solvent_accessibility))

        '''
        elif compute_aromatic_contacts == True:
            isAromatic = rdkit_cage.GetAtomWithIdx(atom_idx_dict[atom_type][dist[1][k]]).GetIsAromatic()
            if isAromatic == True:
                bfactor_info.append(1 / (1 + dist[0][k]))

        '''

        average_cavity_hydrophobicity = np.mean(self.hydrophobicity)
        total_cavity_hydrophobicity = average_cavity_hydrophobicity * self.volume
        logger.info(f"Hydrophobicity method and distance function = {self.hydrophMethod}, {self.distance_function}")
        logger.info(f"Average cavity hydrophobicity = {average_cavity_hydrophobicity:.2f} A^-3")
        logger.info(f"Total cavity hydrophobicity = {total_cavity_hydrophobicity:.2f}")
        return self.hydrophobicity

    def calculate_esp(self, metal_name=None, metal_charge=None, method='eem', max_memory=1e9):
        factor = (8.987551792e9)*(1.602176634e-19)*(1e10) #(Coulomb constant)*(elementary charge)/ Angstrom

        partial_charges = calculate_partial_charges(self.positions, self.atom_names, method=method, metal_name=metal_name, metal_charge=metal_charge)

        if len(self.dummy_atoms_positions) * len(self.positions) * 8 < max_memory:
            dist_matrix = distance_matrix(self.dummy_atoms_positions, self.positions)
            grid_charges = (factor / dist_matrix).dot(partial_charges)
        else:
            grid_charges = np.zeros(len(self.dummy_atoms_positions))
            for idx, dummy_atom in enumerate(self.dummy_atoms_positions):
                dist_matrix = distance_matrix([dummy_atom], self.positions)[0]
                grid_charges[idx] = np.sum(partial_charges * (factor / dist_matrix))

        self.esp_grid = grid_charges

    def calculate_window(self):

        self.window_radius = get_max_escape_sphere(self.positions, self.atom_names)
        return self.window_radius


