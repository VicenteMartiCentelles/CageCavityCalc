from tempfile import mkdtemp
import numpy as np

from CageCavityCalc.input_output import print_to_pdb_file

def calculate_partial_charges_using_ob(positions, atom_names, method='eem'):
    from openbabel import openbabel
    '''
    From Openbabel description:

    eem    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Bultinck B3LYP/6-31G*/MPA
    eem2015ba    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/AIM
    eem2015bm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/MPA
    eem2015bn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/NPA
    eem2015ha    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/AIM
    eem2015hm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/MPA
    eem2015hn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/NPA
    eqeq    Assign EQEq (charge equilibration) partial charges.
    fromfile    Assign charges from file containing {'atom-name', charge} pairs
    gasteiger    Assign Gasteiger-Marsili sigma partial charges
    mmff94       Assign MMFF94 partial charges
    none    Clear all partial charges
    qeq    Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991)
    qtpie    Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007)
    '''
    # openbabel
    try:
        from openbabel import openbabel
    except:
        print("You must install openbabel")

    tmpdir_path = mkdtemp()
    print_to_pdb_file(tmpdir_path + "/temp.pdb", positions, atom_names)

    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("pdb", "mol2")

    # We firstly change pdb -> mol2, because charges are not assigned correctly when used with pdb
    mol = openbabel.OBMol()
    ob_conversion.ReadFile(mol, tmpdir_path + "/temp.pdb")
    ob_conversion.WriteFile(mol, tmpdir_path + "/temp.mol2")
    ob_conversion.SetInAndOutFormats("mol2", "mol2")

    mol = openbabel.OBMol()
    ob_conversion.ReadFile(mol, tmpdir_path + "/temp.mol2")


    ob_charge_model = openbabel.OBChargeModel.FindType(method)

    is_calculated = ob_charge_model.ComputeCharges(mol)

    if is_calculated:
        return ob_charge_model.GetPartialCharges()
    else:
        return None


def calculate_partial_charges(positions, atom_names, method, metal_name=None, metal_charge=None):
    if metal_name is not None:
        list_of_metals = [a for a, atom_name in enumerate(atom_names) if atom_name.title() == metal_name.title()]
        positions = np.delete(positions, list_of_metals, axis=0)
        atom_names = np.delete(atom_names, list_of_metals)
        partial_charges = list(calculate_partial_charges_using_ob(positions, atom_names, method))
        for index in list_of_metals:
            partial_charges.insert(index, metal_charge)
    else:
        partial_charges = list(calculate_partial_charges_using_ob(positions, atom_names, method))
    #print("Sum", np.sum(partial_charges))
    return partial_charges

