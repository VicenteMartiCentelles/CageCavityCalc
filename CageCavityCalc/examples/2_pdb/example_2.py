from CageCavityCalc.CageCavityCalc import cavity

cav = cavity()
cav.read_file("cage_ACIE_2006_45_901.pdb")
cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")
print(cav.volume)
