from main import cavity

cav = cavity()
cav.read_file("cage_1_JACS_2006_128_14120.mol2")
cav.calculate_volume()
cav.calculate_hydrophobicity()
cav.print_to_file("cage_cavity_hydrophobicity.pdb", 'h')
cav.print_to_file("cage_cavity_aromatic_contact.pdb", 'a')
cav.print_to_file("cage_cavity_solvent_accessibility.pdb", 's')
