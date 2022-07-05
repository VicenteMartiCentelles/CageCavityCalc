import numpy as np
import MDAnalysis
from rdkit import Chem

def sum_grid_volume(positions,radius=1.2, volume_grid_size=0.2):
    """
    Calculate the volume of a cavity
    positions: position
    radius:
    volume_grid_size:
    return (float): volume of the cavity
    """
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
    
    
# Assign the hydrophobic values to cage atoms from the library
def assignHydrophobicValuesToCageAtoms(rdkit_cage, hydrophValues, printLevel = 2):
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

    atomTypesValuesListHydrophValues = []
    atomTypesMeanListHydrophValues = []
    atomTypesInfoAtomSymbol = []
    atomTypesInfoAtomGlobalIndex = []
    atomTypesAssignemet = []

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
        else:
            valuesList.append(0)
            print("WARNING: atom ", atomSymbol, i+1, "not found. Assinged 0 as hydropobicity factor.")
        
        atomTypesListHydrophValues.append(valuesList)
        meanValuestList = np.mean(valuesList)
        atomTypesMeanListHydrophValues.append(meanValuestList)
        atomTypesValuesListHydrophValues.append(valuesList)
        atomTypesInfoAtomSymbol.append(atomSymbol)
        atomTypesInfoAtomGlobalIndex.append(i+1)
        atomTypesAssignemet.append(atomTypesList[i])
        if (printLevel == 2): print(atomSymbol, i+1, atomTypesList[i], valuesList, meanValuestList)
    
    return atomTypesMeanListHydrophValues, atomTypesValuesListHydrophValues, atomTypesInfoAtomSymbol, atomTypesInfoAtomGlobalIndex, atomTypesAssignemet
