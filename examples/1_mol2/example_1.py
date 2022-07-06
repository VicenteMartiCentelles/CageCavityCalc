from main import cavity

cav = cavity()
cav.read_file("cage_1_JACS_2006_128_14120.mol2")
cav.calculate_volume()
cav.print_to_file("file.pdb")