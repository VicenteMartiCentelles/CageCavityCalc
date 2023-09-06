
import numpy as np
from rdkit import Chem
from CageCavityCalc.log import logger

# Assign the hydrophobic values to cage atoms from the library
def assignHydrophobicValuesToCageAtoms(rdkit_cage, hydrophValues, printLevel=2):
    listGlobalOut = list()
    listGlobalOut2 = list()
    for i in hydrophValues:
        match = rdkit_cage.GetSubstructMatches(Chem.MolFromSmarts(hydrophValues[i][1]), False, False)
        listOut = list()
        for x in match:
            listOut.append(x[0] + 1)
        listGlobalOut2.append(list(set(listOut)))
        if listOut:
            listGlobalOut.extend(list(set(listOut)))

    numberOfAtoms = rdkit_cage.GetNumAtoms()

    atomTypesList = []
    for i in range(0, numberOfAtoms):
        atomTypesList.append([])

    atomTypesListHydrophValues = []
    for i in range(0, numberOfAtoms):
        atomTypesListHydrophValues.append([])

    atomTypesValuesListHydrophValues = []
    atomTypesMeanListHydrophValues = []
    atomTypesInfoAtomSymbol = []
    atomTypesInfoAtomGlobalIndex = []
    atomTypesAssignemet = []

    for i in range(0, len(listGlobalOut2)):
        if listGlobalOut2[i]:
            for j in listGlobalOut2[i]:
                atomTypesList[j - 1].append(i + 1)

    for i in range(0, numberOfAtoms):
        atomSymbol = rdkit_cage.GetAtomWithIdx(i).GetSymbol()
        valuesList = []
        if len(atomTypesList[i]) > 0:
            valuesList = []
            for j in atomTypesList[i]:
                valuesList.append(hydrophValues[j][2])
        else:
            valuesList.append(0)
            logger.warning(f"WARNING: atom {atomSymbol:}, {i + 1:}, not found. Assinged 0 as hydropobicity factor.")

        atomTypesListHydrophValues.append(valuesList)
        meanValuestList = np.mean(valuesList)
        atomTypesMeanListHydrophValues.append(meanValuestList)
        atomTypesValuesListHydrophValues.append(valuesList)
        atomTypesInfoAtomSymbol.append(atomSymbol.upper()) # we need upper case
        atomTypesInfoAtomGlobalIndex.append(i + 1)
        atomTypesAssignemet.append(atomTypesList[i])
        logger.info(f"{atomSymbol:}, {i + 1:}, {atomTypesList[i]:}, {valuesList:}, {meanValuestList:}")

    return atomTypesMeanListHydrophValues, atomTypesValuesListHydrophValues, atomTypesInfoAtomSymbol, atomTypesInfoAtomGlobalIndex, atomTypesAssignemet



def calc_single_hydrophobicity(distance,Hydroph_Value, distance_function):
    if distance_function == "Audry":
        return Hydroph_Value / (1 + distance)
    elif distance_function == "Fauchere":
        return Hydroph_Value * np.exp(-1 * distance)
    elif distance_function == "Fauchere2":
        return Hydroph_Value * np.exp(-1 / 2 * distance)
    elif distance_function == "OnlyValues":
        return Hydroph_Value