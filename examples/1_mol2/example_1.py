from main import cavity

cav = cavity()
cav.read_file("cage_1_JACS_2006_128_14120.mol2")
cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")
cav.print_to_pymol("cage_cavity.pml")
